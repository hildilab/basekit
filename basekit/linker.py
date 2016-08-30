from __future__ import with_statement
from __future__ import division

import os
import shutil
import itertools
import json
import collections
import xml.etree.ElementTree as ET
import re
import zipfile

from utils import copy_dict, iter_stride, dir_walker
import utils.path
from utils.tool import _, _dir_init, CmdTool, PyTool, ProviMixin
from utils.numpdb import NumPdb, numsele
from utils.mrc import getMrc
import numpy as np

import provi_prep as provi
from spider import LoopCrosscorrel
from pdb import PdbEdit, SplitPdbSSE, LoopDelete, PdbSplit, get_tree
import utils.numpdb as numpdb
import timeit
DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "linker" )


def LINKIT_DIR():
    return os.environ.get("LINKIT_DIR", "")

def LINKIT_DIR2():
    return os.environ.get("LINKIT_DIR2", "")

def LINKIT_CMD():
    return os.path.join( LINKIT_DIR(), "Link_It_dos2n.exe" )

def LINKIT_CMD_mem():
    return os.path.join( LINKIT_DIR2(), "Link_It_dos2n.exe" )

class LinkIt( CmdTool, ProviMixin ):
    """Please upload the PDB file, define the N- and C-terminal stem-residues (e.g.: 10:A, 16:A) and provide the missing sequence in 1-letter code (ACDEF)."""
    args = [
        _( "pdb_file", type="file", ext="pdb", label="PDB File",
            help="The input structure." ),
        _( "res1", type="sele", label="Stem residue 1",
            help="N-terminal stem residue and chain, '123:A'." ),
        _( "res2", type="sele", label="Stem residue 2",
            help="C-terminal stem residue." ),
        _( "seq", type="str", label="Loop sequence",
            help="One-letter code of the linker amino acids." ),
        _( "memdb", type="bool", label="MembraneDB",
            help="Show only results from membrane proteins.", default=False ),
        _( "max_loops", type="int", range=[0, 500], default=100 , step=100,
            advanced=True )
    ]
    out = [
        _( "bin_file", file="{pdb_file.stem}_linker.bin" ),
        _( "txt_file", file="{pdb_file.stem}_linker.txt" ),
        _( "pdb_linker_file", file="{pdb_file.stem}_linker.pdb" ),
        _( "pdb_linker_file2", file="{pdb_file.stem}_linker2.pdb" ),
        _( "pdb_linker_file3", file="{pdb_file.stem}_linker3.pdb" ),
        _( "kos_file", file="{pdb_file.stem}_kos.txt" ),
        _( "json_file", file="{pdb_file.stem}_linker.json" ),
        _( "loop_dir", dir="loops" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it.provi"

    def _init( self, *args, **kwargs ):
        self.res1['chain'] = self.res1['chain'].upper()
        self.res2['chain'] = self.res2['chain'].upper()
        if self.res1['resno'] > self.res2['resno']:
            self.res1, self.res2 = self.res2, self.res1
        if self.memdb:
            self.cmd = [ "wine", LINKIT_CMD_mem(), self.kos_file, self.bin_file, "tp" ]
        else:
            self.cmd = [ "wine", LINKIT_CMD(), self.kos_file, self.bin_file, "tp" ]


    def _pre_exec( self ):
        self._make_kos_file()

    def _post_exec( self ):
        self._fix_linker_pdb( self.pdb_linker_file2 )
        if os.path.getsize(self.pdb_linker_file2) <= 3:
            return
        self._fix_linker_pdb( self.pdb_linker_file3, atoms_only=True )
        self._split_loop_file()
        self._make_linker_json( compact=True )
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            pdb_linker_file3=self.relpath( self.pdb_linker_file3 ),
            json_file=self.relpath( self.json_file )
             )
        
    def _make_kos_file( self ):
        npdb = NumPdb( self.pdb_file, features={
            "phi_psi": False, "sstruc": False, "backbone_only": True
        })

        with open( self.kos_file, "w" ) as fp:
            d = [
                (self.res1, " CA "), (self.res1, " C  "),
                (self.res2, " N  "), (self.res2, " CA ")
            ]
            for sele, atomname in d:
                sele["atomname"] = atomname
                try:
                    coords = npdb.get( 'xyz', **sele )
                    fp.write( "%s\n" % "\n".join(map( str, coords[0] ) ) )
                except Exception as e:
                    self.log(
                        ( "error applying selection with "
                            "residue '%i' and chain '%s'" ) % (
                            sele.get("resno"), sele.get("chain")
                        )
                    )
                    raise e
            fp.write( "%s\n" % self.seq )
            for sele in ( self.res1, self.res2 ):
                fp.write(
                    "%s %s\n" % ( sele.get("chain") or " ", sele["resno"] )
                )

    def _fix_linker_pdb( self, output_file, atoms_only=False, stems=True ):
        backbone = ( ' N  ', ' C  ', ' CA ', ' O  ' )
        chain = self.res1['chain'] or " "
        print chain, self.res1
        with open( self.pdb_linker_file, "r" ) as fp, \
                open( output_file, "w" ) as fp_out:
            for i, line in enumerate( fp ):
                if line.startswith("MODEL"):
                    atom_i = 1
                    if int(line.split () [1])==self.max_loops+1:
                            #fp_out.write( "END" )
                            break
                    else:
                        fp_out.write( line )


                if line.startswith("ENDMDL"):
                    fp_out.write( line )
                if line.startswith("ATOM"):
                    line = line[0:6] + ( "% 5i" % atom_i ) + line[11:]
                    if line[22] == "X":
                        if not stems:
                            continue
                        if line[12:16] not in backbone:
                            continue
                        if line[25] == " ":
                            resnew = int( self.res1['resno'] )
                        else:
                            resnew = int( self.res2['resno'] )
                        resnewp = "%4s" % resnew

                        line = line.replace( "", " " ).replace( "", " " )
                        line = line[0:17] + "GLY" + line[20:]
                        line = line[0:21] + chain + resnewp + " " + line[27:]
                    else:
                        resnew = int( line[22:26] ) + int( self.res1['resno'] )

                        resnewp = "%4s" % resnew
                        line = line = line[0:21] + chain + resnewp + line[26:]
                    atom_i += 1
                    fp_out.write( line )

                    continue
                if not atoms_only:
                    if not line.startswith ("MODEL") and not line.startswith ("ENDMDL"):
                        fp_out.write( line )
            fp_out.write ("END")

    def _split_loop_file( self ):
        PdbSplit(
            self.pdb_linker_file2, output_dir=self.loop_dir,
            backbone_only=True, resno_ignore=[ 1000, 2000 ], zfill=3,max_models=self.max_loops
        )

    def _find_clashes( self ):
        npdb = NumPdb( self.pdb_file, { "backbone_only": True } )

        oripdb=npdb.sele(chain=self.res1['chain'],resno=[self.res1['resno'],self.res2['resno']] ,invert=True)
        noloop= npdb.copy(sele=oripdb)
        nowater=noloop.sele(resname=['HOH','WAT','SOL'],invert=True)
        rest=noloop.copy(sele=nowater)
        rest.write('hurz.pdb')
        


        protein_tree = get_tree(rest['xyz'])
        atom_resno_list = rest.get('resno')

        model_clash_count = {}

        for m, file_path in dir_walker( self.loop_dir, ".*\.pdb" ):
            #oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=True)
            npdb2 = NumPdb( file_path, {"backbone_only": True} )
            loop_tree = get_tree( npdb2['xyz'] )
            k = loop_tree.query_ball_tree( protein_tree, 3 )
            g = [x for x in k if x != []]

            # flatten list of lists
            f = list( itertools.chain(*g) )

            # get unique, sort
            clashatoms = sorted( set( f ) )
            clashes = []

            for i in clashatoms:
                e = atom_resno_list[i]
                if e not in (self.res1['resno']-1, self.res2['resno']+1):
                    clashes.append(e)

            model_no = int(
                utils.path.basename( file_path ).split('_')[0]
            )
            model_clash_count[ model_no ] = len( clashes )


        return model_clash_count
    
         
    def seq_id ( self,seq1,seq2 ):
        si=0
        for sf in range(0, len(seq1),1):
            if seq1[sf]==seq2[sf+1]:
                si+=1
        sqi=(si/len(seq1))*100
        return sqi
    def _make_linker_json( self, compact=False ):
     
        linker_dict = {}
        model_clash_count = self._find_clashes()
        parameter={'res1':self.res1,'res2':self.res2,'sequence':self.seq}
        with open( self.txt_file, "r" ) as fp:
            fp.next()
            fp.next()
            for i, d in enumerate( iter_stride( fp, 4 ), start=1 ):
                if i<=self.max_loops:
                        chain=str(d[3].strip()[-5])
                        stri=re.match(r"([0-9]+)([A-Z]+)", (str(d[3].strip())[-4:]).strip(), re.I)
                        if stri:
                            print stri.group(1)
                            posres1=int(stri.group(1))
                            posres2=int(stri.group(1))+len(self.seq)-1
                        else:
                            posres1=int((str(d[3].strip())[-4:]).strip())
                            posres2=int((str(d[3].strip())[-4:]).strip())+len(self.seq)-1
                        posfield="%s-%s:%s" % (posres1,posres2,chain)
                        linker_dict[ i ] = [
                            float(d[0]), float(d[1]),
                            str(d[2].strip())[1:-1],
                            str(d[3].strip().split() [0]),
                            model_clash_count[ i ],
                            posfield,
                            self.seq_id(self.seq, str(d[2].strip()))     
                        ]
                        
                        
                        
        top_level={
            "params": parameter,
            "linker": linker_dict
        }
        with open( self.json_file, "w" ) as fp:
            if compact:
                
                json.dump(top_level, fp, separators=(',', ':') )
            else:
                json.dump( top_level, fp, indent=4 )

# erstellt eine Statistik zu einem fertig berechneten Datensatz  
# classe SSFEStatistik ( PyTool, ProviMixin ):
#     args = [
#         _("out_datensatz", type="dir")
#     ]
#     
#     def func ( self ) :
        
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
# fuert SSFE fuer einen Datensatz aus 
class SSFEMultiLinkIt( PyTool, ProviMixin ):
    args = [
        _("loop_jobs", type="dir")
    ]   
           
    def func ( self ) :
        
        os.chdir(self.output_dir) 
        
        self.loop_jobs = os.path.abspath(self.loop_jobs)
        
        #enthaelt Liste mit allen Dateien in denen die Proteinen + ihr Loops stehen
        jobDirs = os.listdir (self.loop_jobs)
    
        for jobDir in jobDirs :
            jobDirPath = os.path.join(self.loop_jobs, jobDir)
            inputFiles = os.listdir(jobDirPath)
            jobFile = ''
            pdbFileList = []
           
            for inputFile in inputFiles :
                if  os.path.splitext(inputFile)[1] == '.pdb' :
                    pdbFileList.append(inputFile)
                else :
                    jobFile = inputFile
                    
            # print pdbFileList
            # print jobFile
               
            outJobDir = self.subdir(jobDir + '_Result')
            
            os.chdir(outJobDir)
            
            for pdbFile in pdbFileList : 
                outJobFileDir = self.subdir(os.path.join(outJobDir, os.path.splitext(pdbFile)[0]))
                
                os.chdir(outJobFileDir)
                
                # print '|',pdbFileList, '|',jobFile,'|'
                
                SSFELinkIt(os.path.join(jobDirPath, pdbFile), os.path.join(jobDirPath, jobFile), verbose=self.verbose, debug=True)
        
                os.chdir('..')
                
            os.chdir('..')
        
        #shutil.make_archive(self.output_dir, 'zip' , '.')
        # zipf = zipfile.ZipFile( self.output_dir, 'w', zipfile.ZIP_DEFLATED)
        # for root, dirs, files in os.walk(path) :
        #     for file in files :
        #         ziph.write(os.path.join(root, file))
        # zipf.close()

class SSFELinkIt( PyTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "loop_file", type="file", ext="txt" ),
        _( "extension", type="list", nargs=2, action="append",
           help="int,int", default=[0,0] )
  
    ]
    out = [
        _( "ori_file", file="{pdb_file.stem}_input.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "multi_link_it.provi"
    
    AminoDict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K','ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A','VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    
    def _init( self, *args, **kwargs ):
        
        #falls kein Chain vorhanden in pdb, fuegt 'A' in original_pdb ein
        with open( self.pdb_file, 'r' ) as fp, open( self.ori_file, 'w' ) as fp2:
            for line in fp:
                if line.startswith("ATOM"): #and line[21]==" ":
                    newline = line[:21]+"A"+line[22:]
                    fp2.write(newline)
                else:
                    fp2.write(line)
                 
        #speichert in c den Chain fuer Aufruf       
        c = ''   
        with open( self.ori_file, 'r' ) as fp:
            for line in fp:
                if line.startswith("ATOM") :
                    c = line[21]
                    break
         
        self.multiLinkIts = []
        
        #speichert Loopstart, -ende und Sequenz und Stellen der zu erweiterten Loops
        self.loopTasks = []
        self.loopPosList = []
        #
        loopBrackets = []
        #loopBracketAminos : Stuktur Dict. --> Position:Amino
        loopBracketAminos = {}
        #zaehlt die Uebergebene Loops --> benutzt fuer loopDict
        self.loopCount = 0
        with open( self.loop_file ) as loopTasksInput:
            for line in loopTasksInput:
                if line[0] == '#':
                    continue
                match = re.search("start=(\d*):\s*end=(\d*):\s*(\w*)",line)
                if match:
                    start, end, seq = match.groups()
                    # Wenn Sequenz zu lang, berechne gar nichts
                    if len(seq) > 100 :
                        # self.loopTasks = []
                        # self.loopPosList.append([0, 0])
                        # continue
                        loopBrackets = []
                        loopBracketAminos = {}
                        with open(os.path.join( self.output_dir, "ERROR_LoopSequenzToLong.txt"), 'w') :
                            pass
                        return
                    # fehlerausgabe fuer falsche eingaben
                    if start > end and len(seq) != 0 or start < end and len(seq) == 0 :
                        self.loopTasks = []
                        loopBrackets = []
                        loopBracketAminos = {}
                        with open(os.path.join( self.output_dir, "ERROR_InvalidInput.txt"), 'w') :
                            pass
                        return
                    
                    # loop muss nicht berechnet werden, existiert schon
                    if start > end and len(seq) == 0 :
                        self.loopPosList.append([0, 0])
                        continue
                        
                    self.loopCount += 1        
                
                    #Eingabeliste aller zu bearbeitenden Loops aus Eingabedatei       
                    self.loopTasks.append([int(start)-1 , int(end)+1 , seq])
                    self.loopPosList.append([start, end])
                    
                    #List von Start- und Endwertern der erweiterten Loops
                    startInt = int(start)-1
                    endInt = int(end)+1              
                    loopBrackets.append(([startInt-2,startInt-1,startInt], [endInt, endInt+1, endInt+2]))
                    loopBracketAminos.update(dict.fromkeys(loopBrackets[-1][0]))
                    loopBracketAminos.update(dict.fromkeys(loopBrackets[-1][1]))
        
        self.loopPosList += (7 * [[0, 0]])
        
        #bestimmt Listen mit N und C Terminus der einzelnen Loops
        with open( self.pdb_file, 'r' ) as fp :
            for line in fp :
                if line.startswith("ATOM") :
                    position = int(line[22:26])
                    if position in loopBracketAminos :
                        loopBracketAminos[position] = line[17:20]
        
        #ersetzt 3-Buchstabencode durch 1-Buchstabencode
        for position, aminoTriplet in loopBracketAminos.iteritems() :
            #print position,aminoTriplet
            loopBracketAminos[position] = SSFELinkIt.AminoDict[aminoTriplet]
            
        #print loopBracketAminos                      
        
        for i in range (4) :
            loopTasksList = []
            for start,end,seq in self.loopTasks :
                for j in range(i) :
                    seq = loopBracketAminos[start - j ] + seq + loopBracketAminos[end + j]    
                loopTasksList.append([str(start - i) + ':' + c, str(end + i) + ':' + c, seq])
             
            #print loopTasksList
               
            self.multiLinkIts.append(MultiLinkIt(self.ori_file, input = loopTasksList,
                    **copy_dict( kwargs, run=False,
                                output_dir=self.subdir("link_it_%i,%i" % (i,i) ), verbose=True, debug=True )))
            
                     
    def func( self ):
        if len(self.multiLinkIts) == 0 :
            return
        
        for m in self.multiLinkIts : 
            m()
        
        #loopDict enthaelt fuer jeden Loop die 3 besten berechneten Loops pro Erweiterung:(0,0);(1,1);(2,2);(3,3)
        loopDict = {}
        #zum Laden der einzelnen Loop-json zur Auswahl der 3 besten Loops nach Score und Anzahl Clashes
        singleLoopDict = {}
        
        # print self.loopCount
        
        #erstellt Dict mit den 3 besten berechneten Loops pro Loop 
        for i in range(4) :
            for j in range(self.loopCount) :
                jsonPfad = os.path.join( self.output_dir,"link_it_%i,%i/link_it_%i" % (i,i,j) )+ "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker.json"
                if os.path.isfile(jsonPfad) :
                    with open(os.path.join( self.output_dir,"link_it_%i,%i/link_it_%i" % (i,i,j) )+ "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker.json", 'r') as fp :
                        singleLoopDict = json.load(fp)
                        loop = loopDict.get(j, [])
                        #print type(loopDict)
                        #print (os.path.join( self.output_dir,"link_it_%i,%i/link_it_%i" % (i,i,j) )+ "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker.json")
                    for resultIndex in range(1, 4) :
                        result = singleLoopDict["linker"][str(resultIndex)]
                        loop.append(["%i,%i;%i;%i" % (i,i,j,resultIndex), result[1], result[4], result[2], result[3], result[5]])
                        loopDict[j] = loop
                        
        #Ergebnis wird in BestResultsDict.json ausgegeben                
        with open(os.path.join(self.output_dir, "BestResultsDict.json"), 'w') as fp :
            json.dump(loopDict, fp)
        # print loopDict
        
        #topTenResult entaelt die Loopinformationen der 10 besten gesamt Ergebnisse
        # --> benutzt fuer zusammenbauen der 10 besten pdb-files
        self.topTenResult = {}
        
        #sortiert nach hoestem Score die Liste der 3 besten berechneten Loops pro Loop (12 Ergebnisse pro Loop) und dreht sie um, d.h bester Score zu erst
        for key in loopDict :
            loopDict[key].sort(key=lambda loop: (loop[1], loop[2]))
            loopDict[key].reverse()
        
        #erzeugt TopTen Liste von 1-10
        minNumResults = min(10, min(map(len, loopDict.values())))
        for i in range(1,minNumResults + 1) :
            resultList = []
            for key, loopList in loopDict.iteritems() :
                resultList.append(loopList[i-1])
             
            self.topTenResult[i] = resultList
            
        #Ergebnis wird in TopTenResult.json ausgegeben
        with open(os.path.join(self.output_dir, "TopTenResult.json"), 'w') as fp :
            json.dump(self.topTenResult, fp)
            
        self.insertLoop()
                   
                                          
    #fuegt Loops in pdb-file ein
    #--> Ergebnis sind die 10 besten pdb-files
    def insertLoop( self ):
        
        #shutil.copyfile('../../../index.html', os.path.join(self.output_dir, "index.html" ))
        
        templateHtmlString = ''
        
        for key in self.topTenResult :
            loopLineList = []
            for loop in self.topTenResult[key] :
                extension, loopName, numLoop = loop[0].split(';')
                
                loopName = int(loopName)
                numLoop = int(numLoop)
                #da extension ein Doupel ist (z.B. 0,0), wird es aufgeteilt 
                extensionSingle1, extensionSingle2 = extension.split(',')
                extensionSingle1 = int(extensionSingle1)
                extensionSingle2 = int(extensionSingle2)
                
                loopFile = []
                with open(self.subdir("link_it_%i,%i/link_it_%i" % (extensionSingle1,extensionSingle2,loopName)) + "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker2.pdb", 'r') as fp :
                    loopFile = fp.readlines()
                
                modelBegin = -1
                modelEnd = -1
                for i, line in enumerate(loopFile) :
                    if line[0:5] == 'MODEL' and line.split()[1] == str(numLoop) :
                        modelBegin = i
                    if line[0:6] == 'ENDMDL' and modelBegin != -1 :
                        modelEnd = i
                        break
                        
                loopLineList.append(loopFile[modelBegin + 5 : modelEnd - 5])
                
                #print modelBegin
                #print modelEnd
                #print loopLineList 
        
                       
            new_pdb_file = (os.path.splitext(self.ori_file)[0] + '_out_%i.pdb' %key)
            with open( self.ori_file, 'r') as fp, open (new_pdb_file, 'w') as fp2 :
                for line in fp :
                    #Filter ANISOU-Zeilen raus, falls vorhanden und kopiert alle anderen Zeilen
                    if line[0:6] != 'ANISOU' :
                        fp2.write(line)
            for loop in loopLineList :
                    
                startPos = int(loop[0].split()[5])
                endPos = int(loop[-1].split()[5])
                # print startPos
                # print endPos
                # print loop

                #Fuegt Loop an richtiger Stelle in Datei ein
                rawPdbFile = []
                with open( new_pdb_file, 'r') as fp :
                    rawPdbFile = fp.read().splitlines(True)
                with open (new_pdb_file + '.tmp', 'w') as fp2 :
                        
                    for i, line in enumerate(rawPdbFile) :
                        #print i, line[0:4], line[22:26]
                        if ((line[0:3]) == 'TER' and int(line[22:26]) == startPos -1 ) or \
                           ((line[0:4]) == 'ATOM' and (int(line[22:26]) in range(startPos,endPos+1))) or \
                           ((line[0:3]) == 'TER' and (int(line[22:26]) in range(startPos,endPos+1))) :
                            for loopLine in loop :
                                fp2.write(loopLine)
                            loop = []
                        elif (((line[0:4]) == 'ATOM' and int(line[22:26]) == endPos +1) and \
                             ((rawPdbFile[i-1][0:4]) == 'ATOM' and int(rawPdbFile[i-1][22:26]) == startPos -1)) :
                            for loopLine in loop :
                                fp2.write(loopLine)
                            loop = []
                            fp2.write(line)
                        else :
                            fp2.write(line)
                
                shutil.move(new_pdb_file + '.tmp', new_pdb_file)
                    
                    
                #Neu-Nummerierung von out_pdb             
                i = 0
                with open( new_pdb_file, 'r') as fp, open (new_pdb_file + '.tmp', 'w') as fp2 :
                    for line in fp:
                        if line.startswith("ATOM"):
                            i= i+1
                            newline = line[:6] + " "*(5-len(str(i))) + str(i) + line[11:]
                            fp2.write(newline)
                        else:
                            fp2.write(line)
                            
                shutil.move(new_pdb_file + '.tmp', new_pdb_file)
            
            
            
            #html.index erstellen            
            shutil.copyfile('../../../ngl.embedded.min.js', os.path.join(self.output_dir, "ngl.embedded.min.js" ))
            # tmpl_name="ngl.embedded.min.js"
            # tmpl_file = os.path.join( tmpl_dir, tmpl_name )
            templateHtmlString += '"' + os.path.basename(new_pdb_file) + '",\n'
            # print templateHtmlString
            
            templateHtml = ''
            with open ('../../../index.html', 'r') as fp :
                templateHtml = fp.read()
            
            #print self.loopPosList
            templateHtml = templateHtml.format(namepdbfiles = templateHtmlString, loopPos1 = '"'+ str(self.loopPosList[0][0]) + '-' + str(self.loopPosList[0][1])+'"',
                                               loopPos2 = '"'+ str(self.loopPosList[1][0]) + '-' + str(self.loopPosList[1][1])+'"',
                                               loopPos3 = '"'+ str(self.loopPosList[2][0]) + '-' + str(self.loopPosList[2][1])+'"',
                                               loopPos4 = '"'+ str(self.loopPosList[3][0]) + '-' + str(self.loopPosList[3][1])+'"',
                                               loopPos5 = '"'+ str(self.loopPosList[4][0]) + '-' + str(self.loopPosList[4][1])+'"',
                                               loopPos6 = '"'+ str(self.loopPosList[5][0]) + '-' + str(self.loopPosList[5][1])+'"',
                                               loopPos7 = '"'+ str(self.loopPosList[6][0]) + '-' + str(self.loopPosList[6][1])+'"')
            # print templateHtml
            
            
            with open (os.path.join(self.output_dir, "index.html" ), 'w') as fp :
                fp.write(templateHtml)
   
                      
class MultiLinkIt( PyTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "input", type="list", nargs=3, action="append",
           help="sele,sele,str", default=None ),
        _( "names", type="list", nargs="*", default=None )
    ]
    out = [
        _( "ori_file", file="{pdb_file.stem}_out.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "multi_link_it.provi"

    def _init( self, *args, **kwargs ):
        shutil.copyfile(self.pdb_file, self.ori_file)
        self.pdb_file = self.ori_file
        self.link_it_list = []
        for i, linker_args in enumerate( self.input ):
            #print linker_args
            res1, res2, seq = linker_args
            #print 'res1:  ' , res1
            #print 'res2:  ' , res2
            #print 'seq:  ' , seq
            link_it = LinkIt(
                self.pdb_file, res1, res2, seq,
                **copy_dict( kwargs, run=False,
                             output_dir=self.subdir("link_it_%i" % i) )
            )
            self.output_files += link_it.output_files
            self.sub_tool_list.append( link_it )
            self.link_it_list.append( link_it )
            

    def func( self ):
        self.log( "%i linkit runs" % len( self.input ) ) 
        for i, link_it in enumerate ( self.link_it_list ):
            link_it()
            
            #print link_it.res1["resno" ] +1
            #print link_it.res2["resno" ] -1
            
            #print self.subdir("link_it_%i" % i + "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_linker2.pdb")
            #self.insertLoop(i)
            
            
    def _post_exec( self ):
        linker_list = []
        for i, link_it in enumerate( self.link_it_list ):
            linker_list.append({
                "i": i,
                "name": self.names[i] if self.names else "Linker",
                "json_file": self.relpath( link_it.json_file ),
                "pdb_linker_file3": self.relpath( link_it.pdb_linker_file3 )
            })
        
        # self._make_provi_file(
        #     use_jinja2=True,
        #     pdb_file=self.relpath( self.pdb_file ),
        #     linker_list=linker_list
        # )

class LinkItDensity( PyTool ):
    """Please denote the resolution of your map. Then define the stem-residues (e.g.: 10:A, 16:A) and provide the missing sequence in 1-letter code (ACDEF). """
    args = [
        _( "pdb_file", type="file", ext="pdb", label="PDB File",
            help="The input structure." ),
        _( "mrc_file", type="file", ext="mrc", label="MRC File",
            help="The input density." ),
        _( "res1", type="sele", label="Stem residue 1",
            help="N-terminal stem residue and chain, '123:A'." ),
        _( "res2", type="sele", label="Stem residue 2",
            help="C-terminal stem residue." ),
        _( "seq", type="str", label="Sequence",
            help="One-letter code of the linker amino acids." ),
        _( "resolution", type="float", fixed=True , range=[0.5,20], default=5,
            precision=1, label="Map resolution",
            help="Used for filtering the map." ),
        _( "max_loops", type="int", range=[0, 500], default=100 , step=100,
            advanced=True )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    #start=timeit.timeit()
    def _init( self, *args, **kwargs ):
        if self.res1['resno'] > self.res2['resno']:
            self.res1, self.res2 = self.res2, self.res1
        self.link_it = LinkIt(
            self.edited_pdb_file, self.res1, self.res2, self.seq,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("link_it") )
        )
        
        self.loop_correl = LoopCrosscorrel(
            self.mrc_file, self.pdb_file,
            self.link_it.pdb_linker_file2,
            self.link_it.txt_file,
            self.res1, self.res2, len(self.seq),
            self.resolution,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("loop_correl"),
                max_loops=self.max_loops,
            )
        )
        self.output_files += list( itertools.chain(
            self.link_it.output_files,
            self.loop_correl.output_files,
        ))
        self.sub_tool_list.extend( [
            self.link_it, self.loop_correl
        ] )
        

    def func( self ):
        boxsizex = getMrc(self.mrc_file, 'nx' )
        boxsizey=getMrc(self.mrc_file,'ny' )
        boxsizez=getMrc(self.mrc_file,'nz' )
        originx = getMrc(self.mrc_file, 'nxstart' )*-1
        originy = getMrc(self.mrc_file, 'nystart' )*-1
        originz = getMrc(self.mrc_file, 'nzstart' )*-1
        xorg=getMrc(self.mrc_file,'xorg')*-1
        yorg=getMrc(self.mrc_file,'yorg')*-1
        zorg=getMrc(self.mrc_file,'zorg')*-1
        size = getMrc(self.mrc_file, 'xlen' )
        pixelsize = size / boxsizex
        if originx!=0 and xorg==0 :
            shx = (originx -(boxsizex/2)) * pixelsize
            shy = (originy -(boxsizey/2)) * pixelsize
            shz = (originz -(boxsizez/2)) * pixelsize
        else:
            shx = ((xorg  / pixelsize)-(boxsizex/2))*pixelsize
            shy = ((yorg/ pixelsize)-(boxsizey/2))*pixelsize
            shz = ((zorg / pixelsize)-(boxsizez/2))*pixelsize
        PdbEdit(
            self.pdb_file, shift=[ shx, shy, shz ]
        )

        self.link_it()
        self.loop_correl()
        #print len(self.seq)
    def _post_exec( self ):
        self._make_correl_json()
        #self._make_fixed_linker()

    def _make_correl_json( self, compact=False ):
        li = self.link_it.json_file
        cc = self.loop_correl.spider_crosscorrelation.crosscorrel_json
        with open( li, "r" ) as fp:
            li_dict = json.load(
                fp, object_pairs_hook=collections.OrderedDict
            )
        with open( cc, "r" ) as fp:
            cc_dict = json.load(
                fp, object_pairs_hook=collections.OrderedDict
            )
        li_dict_li={}   
        li_dict_li=li_dict["linker"]
        linker_correl_dict = {}
        parameter={'res1':self.res1,'res2':self.res2,'sequence':self.seq,'resolution':self.resolution}
        for k, v in cc_dict.iteritems():
            linker_correl_dict[k] = [ v ] + li_dict_li[k]
        top_level={
            "params": parameter,
            "linker": linker_correl_dict
        }
        with open( self.linker_correl_file, "w" ) as fp:
            if compact:
                json.dump( top_level, fp, separators=(',', ':') )
            else:
                json.dump( top_level, fp, indent=4 )

    #end = timeit.timeit()
    #print "zEdIT", end - start

class LnkItVali(PyTool, ProviMixin):
    args = [
        _( "dataset_dir", type="dir" ),
        _( "max_loops", type="int", range=[0, 200], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    def func( self ):
        for fn in os.listdir(self.dataset_dir):
            for files in os.listdir(self.dataset_dir+"/"+fn):
                if files.endswith(".pdb"):
                    SplitPdbSSE(self.dataset_dir+"/"+fn+"/"+files, output_dir=self.dataset_dir+"/"+fn )
        for fn in os.listdir(self.dataset_dir):
            print fn
            for files in os.listdir(self.dataset_dir+"/"+fn+"/pieces"):
                if files.endswith(".json"):
                    jason= os.path.join(self.dataset_dir,  "%s/%s/%s" % (fn,'pieces',files))
                    pdb=os.path.join(self.dataset_dir,"%s/%s.%s" % (fn, files[0:4],'pdb'))
                    emmap=self.dataset_dir+"/"+fn+"/map/"+fn[0:3].lower()+"_"+fn[4:8]+".map"
                    emmap=os.path.join(self.dataset_dir,"%s/%s/%s_%s.%s" % (fn,'map',fn[0:3].lower(),fn[4:8],'map'))
                    xml_file=os.path.join(self.dataset_dir, "%s/%s/%s.%s" % (fn,'header',fn.lower(),'xml'))
                    print pdb
                    print emmap
                    print xml_file
                    with open(jason,"r") as hz:
                        js=json.load(hz)
                        print js["3"]
                        res1=js["1"]
                        res2=js["2"]
                        seq=js["3"]
                        chain=js["4"]
                        stem1=str(res1)+":"+chain
                        stem2=str(res2)+":"+chain
                        del1=res1+1
                        del2=res2-1
                        print stem1
                    with open (xml_file, "r") as rs:
                        tree = ET.parse(rs)
                        root = tree.getroot()
                        resu=float(root.find("processing/reconstruction/resolutionByAuthor").text)
                        print resu
                    LoopDelete (pdb,chain,del1,del2,output_dir=self.dataset_dir+"/"+fn+"/pieces/"+files[:-5]   )

                    if (len(seq)>2) and(len(seq)<20):
                        try:
                            LinkItDensity (
                                self.dataset_dir+"/"+fn+"/pieces/"+files[:-5]+"/noloop.pdb", emmap,stem1 ,stem2, seq ,resu,  output_dir=self.dataset_dir+"/"+fn+"/pieces/"+files[:-5]
                            )
                        except:
                            print "linikt error"
                            continue
                    else:
                        print "length"
                        continue
#   Analyse
#found the original fragment?
class CutPDB (PyTool):
    args = [
    _( "pdb_file", type="file"),
    _( "mrc_file", type="file"),
    _( "resolution", type="float", range=[1, 10], step=0.1 ),
        ]
    def func( self, *args, **kwargs ):

        npdb=numpdb.NumPdb( self.pdb_file)

        for z, numa in enumerate (npdb.iter_chain()):
            cha=numa.get('chain')[0]
            rel=numa.get('resno')[-1]
            #cha='P'
            if cha=='P':
                print rel
                print cha
                os.mkdir(cha)
                for x in range (18,35,1):
                    #ch=numa.get('chain')[0]
                    gum= "%s/%s" % (cha,x)
                    os.mkdir(gum)
                    for i,numa in enumerate (npdb.iter_resno(chain=cha)):



                        res1=numa.get('resno')[0]
                        res2=res1+x-1
                        if res1 % 7 == 0:
                            print res1
                            test=npdb.copy(resno=[res1,res2],chain=cha)
                            oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=True)
                            if rel-x+1>=res1:
                                di="%s/%s_%s" % (gum,res1,res2)
                                os.mkdir(di)
                                name= "%s/%s_%s.%s" % (di,res1,res2,'pdb')
                                name2= "%s/%s_%s_%s.%s" % (di,'ganz',res1,res2,'pdb')
                                print name2
                                test.write(name)
                                npdb.write(name2,sele=oripdb)
                                seq=test.sequence()
                                ires1="%s:%s" % (res1-1,cha)
                                ires2="%s:%s" % (res2+1,cha)
                                print ires1, ires2
                                print 'sequence',seq
                                #print 'richtige?',numa.sequence()
                                try:
                                    LinkItDensity (name2,self.mrc_file,ires1,ires2,seq,self.resolution, output_dir=di)
                                except:
                                    print 'linikt error'
                                #    neuen loop suchen
                else:
                    continue
class CutPDB2 (PyTool):
    args = [
    _( "pdb_file", type="file"),
    _( "mrc_file", type="file"),
    _( "resolution", type="float", range=[1, 10], step=0.1 )
        ]
    def func( self, *args, **kwargs ):

        npdb=numpdb.NumPdb( self.pdb_file)

        for z, numa in enumerate (npdb.iter_chain()):
            cha=numa.get('chain')[0]
            rel=numa.get('resno')[-1]
            #cha='P'
            if cha=='A':
                print rel
                print cha
                os.mkdir(cha)
                for x in range (5,35,1):
                    #ch=numa.get('chain')[0]
                    gum= "%s/%s" % (cha,x)
                    os.mkdir(gum)
                    for i,numa in enumerate (npdb.iter_resno(chain=cha)):



                        res1=numa.get('resno')[0]
                        res2=res1+x-1
                        if res1 % 3 == 0:
                            print res1
                            test=npdb.copy(resno=[res1,res2],chain=cha)
                            oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=True)
                            if rel-x+1>=res1:
                                di="%s/%s_%s" % (gum,res1,res2)
                                os.mkdir(di)
                                name= "%s/%s_%s.%s" % (di,res1,res2,'pdb')
                                name2= "%s/%s_%s_%s.%s" % (di,'ganz',res1,res2,'pdb')
                                print name2
                                test.write(name)
                                npdb.write(name2,sele=oripdb)
                                seq=test.sequence()
                                ires1="%s:%s" % (res1-1,cha)
                                ires2="%s:%s" % (res2+1,cha)
                                print ires1, ires2
                                print 'sequence',seq
                                #print 'richtige?',numa.sequence()
                                try:
                                    LinkItDensity (name2,self.mrc_file,ires1,ires2,seq,self.resolution, output_dir=di)
                                except:
                                    print 'linikt error'
                                #    neuen loop suchen
                else:
                    continue

class AnalyseLiniktRun( PyTool , ProviMixin):
    args = [
    _( "dataset_dir", type="dir" )
    ]
    out = [
     _( "linker_correl_file", file="linker_correl.json" )
    ]
    def func (self):
        for fn in os.listdir(self.dataset_dir):
            #print fn
            for pn in os.listdir(self.dataset_dir+"/"+fn+"/pieces"):

                if os.path.isdir(os.path.join(self.dataset_dir,fn,'pieces',pn)):
                    try:
                        #print pn
                        pf=os.path.join(self.dataset_dir,fn,'pieces',pn,"loop_correl","crosscorrelation")
                        ps=os.path.join(self.dataset_dir,fn,'pieces',pn,"link_it")
                        lj= "%s/%s.json" % (ps,"edited_linker")
                        pdb=pn[0:4].lower()
                        #print pdb
                        with open(lj,"r") as hz:
                            js=json.load(hz)
                            #print js["1"][3]
                            #print "hallo"
                            for item in js:
                                if js[item][3]== pdb:
                                    print "YEEHA"
                                else:
                                    #print "ooh"
                                    continue
                    except:
                        continue
                else:
                    continue
                    #print "%s:%s" % (pn, 'is no directory')




#/home/jochen/work/fragfit/validation/dataset/EMD-5249/pieces/3IZM_B_705-708/loop_correl/crosscorrelation
