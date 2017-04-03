from __future__ import with_statement
from __future__ import division

import os
import shutil
import itertools
import json
import collections
import matplotlib.pyplot as pyplot
from pylab import rcParams
rcParams['figure.figsize'] = 16,9
import xml.etree.ElementTree as ET
import re
import gzip
import zipfile
import time

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
index_html_file = os.path.join( TMPL_DIR, "index.html" )
ngl_js_file = os.path.join( TMPL_DIR, "ngl.js" )
downloader_js_file = os.path.join( TMPL_DIR, "downloader.js" )

def LINKIT_DIR():
    return os.environ.get("LINKIT_DIR", "")

def LINKIT_DIR2():
    return os.environ.get("LINKIT_DIR2", "")

def LINKIT_DIR3():
    return os.environ.get("LINKIT_DIR3", "")

def LINKIT_CMD():
    return os.path.join( LINKIT_DIR(), "Link_It_dos2n.exe" )

def LINKIT_CMD_mem():
    return os.path.join( LINKIT_DIR2(), "Link_It_dos2n.exe" )

def LINKIT_CMD_GPCR():
    return os.path.join( LINKIT_DIR3(), "Link_It_dos2n.exe" )

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
        _( "GPCRdb", type="bool", label="GPCRDB",
            help="Show only results from GPCR proteins.", default=False ),
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
        elif self.GPCRdb:
            self.cmd = [ "wine", LINKIT_CMD_GPCR(), self.kos_file, self.bin_file, "tp" ]
        else:
            self.cmd = [ "wine", LINKIT_CMD(), self.kos_file, self.bin_file, "tp" ]


    def _pre_exec( self ):
        self._make_kos_file()

    def _post_exec( self ):
        if os.path.isfile(self.pdb_linker_file):
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
        else:
            with open(os.path.join( self.output_dir, "ERROR_no_loops.txt"), 'w') :
                pass

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
        #print chain, self.res1
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
                    line = line[0:6] + ( "% 5i" % atom_i ) + line[11:]#60]+ " "+ line[61:]
                    if line[60]!=" ":
                        line = line[0:60]+line[61:]
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
        #print 'jkh',rest['xyz']
        atom_resno_list = rest.get('resno')

        model_clash_count = {}

        for m, file_path in dir_walker( self.loop_dir, ".*\.pdb" ):
            #oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=True)
            #print m, file_path
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
                            #print stri.group(1)
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

#erstellt eine Statistik zu einem fertig berechneten Datensatz
class SSFEStatistic ( PyTool, ProviMixin ):
    args = [
        _("out_dataSet", type="dir"),
        _( "GPCRscore", type="int", default=20 ),
        _( "Speciesscore", type="int", default=1000 )
    ]

    def func ( self ) :

        statTopTenResultList = []
        statBestResultsList = []
        statAvgScoreList = []
        statResultLoopFoundList = []
        statSortedPdbLoopDictList = []

        filename = 'RelResultLoopsList.json'
        resultLoopFoundList = []

        #Ergebnisse fuer alle berechneten Spezies zusammenfassen
        os.chdir(self.out_dataSet)
        for result_dir in [d for d in os.listdir(self.out_dataSet) if os.path.isdir(d)] :
            for template_dir in [d for d in os.listdir(result_dir) if os.path.isdir(os.path.join(result_dir,d))]:
                path = os.path.join(result_dir, template_dir)
                if filename in os.listdir(path) :
                    #print filename
                    with open (os.path.join(path, filename), 'r') as fp :
                        resultLoopList = json.load(fp)
                    resultLoopList.insert(0,result_dir)
                    resultLoopList.insert(1,template_dir)
                    resultLoopFoundList.append(resultLoopList)
                    #print resultLoopFoundList

        with open(os.path.join(self.output_dir, "ResultLoopFoundList.json"), 'w') as fp :
            json.dump(resultLoopFoundList, fp)

        with open (os.path.join(self.output_dir, 'ResultLoopFoundList.csv'), 'w' ) as fp7 :
            fp7.write( 'Timestamp, Template, Name, ICL1, ICL2, ICL3, ECL1, ECL2, ECL3, TM7_H8\n')
            for result in resultLoopFoundList :
                fp7.write(result[0] + ',' + result[1] + ',' + result[2] + ',' + str(result[3][1]) + ',' + str(result[4][1]) + ',' + str(result[5][1]) + ',' + str(result[6][1]) + ','
                          + str(result[7][1]) + ',' + str(result[8][1]) + ',' + str(result[9][1]) + '\n')


        for subdir, dirs, files in os.walk(self.out_dataSet):
            if "BestResultsDict.json" in files and "TopTenResult.json" in files and "AvgScore.json" and "SortedPdBLoopDict.json" in files :
                with open (os.path.join(subdir, "BestResultsDict.json"), 'r') as fp :
                    bestResultsList = json.load(fp)
                    statBestResultsList += [x for sublist in bestResultsList.values() for x in sublist]
                with open (os.path.join(subdir, "TopTenResult.json"), 'r') as fp2 :
                    topTenResultList = json.load(fp2)
                    statTopTenResultList += [x for sublist in topTenResultList.values() for x in sublist]
                with open (os.path.join(subdir, "AvgScore.json"), 'r') as fp3 :
                    avgScoreList = json.load(fp3)
                    statAvgScoreList += [x for x in avgScoreList]
                with open (os.path.join(subdir, "SortedPdBLoopDict.json"), 'r') as fp5 :
                    sortedPdbLoopDictList = json.load(fp5)
                    statSortedPdbLoopDictList += [x for sublist in sortedPdbLoopDictList.values() for x in sublist]



        with open (os.path.join(self.output_dir, "ResultLoopFoundList.json"), 'r') as fp4 :
            resultLoopFoundList = json.load(fp4)
            statResultLoopFoundList += [x for x in resultLoopFoundList]

        #print statResultLoopFoundList

        countLoopFoundList = {'ICL1':[0,0], 'ICL2':[0,0], 'ICL3':[0,0], 'ECL1':[0,0], 'ECL2':[0,0], 'ECL3':[0,0], 'TM7_H8':[0,0]}
        for loop in statResultLoopFoundList :
            #print loop
            for res in loop[3:] :
                #print res
                if res[1] :
                    countLoopFoundList[res[0]][0] +=1
                else :
                    countLoopFoundList[res[0]][1] +=1

        #print countLoopFoundList

        countLoopSpeziesFoundList = {key : [set(), set()] for key in ['ICL1', 'ICL2', 'ICL3', 'ECL1', 'ECL2', 'ECL3', 'TM7_H8']}
        for loop in statResultLoopFoundList :
            for res in loop[3:] :
                if res[1] :
                    countLoopSpeziesFoundList[res[0]][0].add(loop[2])
                else :
                    countLoopSpeziesFoundList[res[0]][1].add(loop[2])

        countSpezies = set()

        for loop in statResultLoopFoundList :
            countSpezies.add(loop[2])

        countSpeziesNum = len(countSpezies)

        for key in countLoopSpeziesFoundList :
            countLoopSpeziesFoundList[key] = [len(countLoopSpeziesFoundList[key][0]), len(countLoopSpeziesFoundList[key][1])]



        with open (os.path.join(self.output_dir, 'CountLoopResults.csv'), 'w' ) as fp6 :
            fp6.write('Loop, found/total, found/species, not found/total, not found/species, could be found, species found\n' )
            for key, result in countLoopFoundList.iteritems():
                fp6.write(key + ',' + str(result[0]) + ',' + str(countLoopSpeziesFoundList[key][0]) + ',' + str(result[1]) + ',' + str(countLoopSpeziesFoundList[key][1]) + ','
                          + str(result[0] + result[1]) + ',' + str(countSpeziesNum) + '\n')

        # print countLoopFoundList
        # print '================================'
        # print countLoopSpeziesFoundList

        #print statTopTenResultList
        # print '----------------'
        # print statSortedPdbLoopDictList

        for x in statTopTenResultList :
            extension, loopName, whichLoop, datenbank = x[0].split(';')
            exStart, exEnd = extension.split(',')
            x[:] = [int(exStart), int(exEnd), int(loopName), int(whichLoop)] + x[1:]
            loopLen = len(x[6])
            x.append(loopLen)
            oriLoopLen = x[13] - x[0] - x[1]
            x.append(oriLoopLen)
            x.append(int(datenbank))

        for x in statBestResultsList :
            extension, loopName, whichLoop, datenbank = x[0].split(';')
            exStart, exEnd = extension.split(',')
            x[:] = [int(exStart), int(exEnd), int(loopName), int(whichLoop)] + x[1:]
            loopLen = len(x[6])
            x.append(loopLen)
            oriLoopLen = x[13] - x[0] - x[1]
            x.append(oriLoopLen)
            x.append(int(datenbank))

        print statSortedPdbLoopDictList

        for x in statSortedPdbLoopDictList :
            extension, loopName, whichLoop, datenbank = x[0].split(';')
            exStart, exEnd = extension.split(',')
            x[:] = [int(exStart), int(exEnd), int(loopName), int(whichLoop)] + x[1:]
            loopLen = len(x[6])
            x.append(loopLen)
            oriLoopLen = x[13] - x[0] - x[1]
            x.append(oriLoopLen)
            x.append(int(datenbank))
            newscore = x[4]
            # print newscore
            if x[11]:#Spezies gefunden
                newscore = newscore/self.Speciesscore
            elif x[9]:#GPCR gefunden
                newscore = newscore/self.GPCRscore
            x[4] = newscore
            # print'-------'
            # print x[4]
            # print 'next'

        # print '-------------------'
        print statSortedPdbLoopDictList
        # [0][1] Erweiterung, [2] Loop, [3] Loopauswahl, [4] Score, [5] Clashes, [6] Sequenz, [7] Template, [8] Position(Datenbank), [9] GPCR True/False, [10] SequenzIdentitaet, [11] Spezies True/False, [12] memDB True/False, [13] Looplaenge, [14] originalLooplaenge, [15] datenbank 0=GPCRDB 1=memDB


        with open (os.path.join(self.output_dir, "StatBestResultsList.json"),'w' ) as fp :
            json.dump(statBestResultsList, fp)

        with open (os.path.join(self.output_dir, "StatTopTenResultList.json"), 'w') as fp2 :
            json.dump(statTopTenResultList, fp2)

        with open (os.path.join(self.output_dir, "StatAvgScoreList.json"), 'w') as fp3 :
            json.dump(statAvgScoreList, fp3)

        with open (os.path.join(self.output_dir, "StatResultLoopFoundList.json"), 'w') as fp4 :
            json.dump(statResultLoopFoundList, fp4)

        with open (os.path.join(self.output_dir, "StatSortedPdbLoopDictList.json"), 'w') as fp5 :
            json.dump(statSortedPdbLoopDictList, fp5)


        stat1Dir = self.subdir('AvgScoreVsLoopLengthPerLength')
        stat2Dir = self.subdir('RateVsExtension')
        stat3Dir = self.subdir('ScoreVsLoopLengthPerExtension')
        stat4Dir = self.subdir('ClashVsLoopLengthPerExtension')
        #stat5Dir = self.subdir('AvgScoreVsLoopLenghtPerModel')
        stat6Dir = self.subdir('DataspaceAndGPCR')

        #wie oft wurde welche Datenbank genutz
        gpcrDB = 0
        memDB =0

        for j, loop in enumerate(statSortedPdbLoopDictList) :
            if loop[15] == 0:
                gpcrDB +=1
            elif loop[15] == 1 :
                memDB +=1

        self.dataspace = [0,1]
        self.numDataspace = [gpcrDB, memDB]

        pyplot.close()
        pyplot.xlabel('dataspace')
        pyplot.ylabel('rate')
        pyplot.xticks(range(2))
        pyplot.title('rate per dataspace')
        pyplot.bar(self.dataspace, self.numDataspace, color = '#BDBDBD')
        pyplot.savefig(os.path.join(stat6Dir, 'RateDataspace.svg'), format='svg')

        # Verhaltniss von gefundener Spezies und GPCR
        FalseSpeziesFalseGpcr = 0
        TrueSpeziesTrueGpcr = 0
        FalseSpeziesTrueGpcr = 0
        TrueSpeziesFalseGpcr = 0

        for j, loop in enumerate(statSortedPdbLoopDictList) :
            if loop[11] == False and loop[9] == False :
                FalseSpeziesFalseGpcr +=1
            elif loop[11] == True and loop[9] == True :
                TrueSpeziesTrueGpcr +=1
            elif loop[11] == False and loop[9] == True :
                FalseSpeziesTrueGpcr +=1
            elif loop[11] == True and loop[9] == False :
                TrueSpeziesFalseGpcr +=1

        self.speziesVsGpcr = [0,1,2,3]
        self.numSpeziesVsGpcr = [FalseSpeziesFalseGpcr, TrueSpeziesTrueGpcr, FalseSpeziesTrueGpcr, TrueSpeziesFalseGpcr]

        print 'SpeziesAndGPCR:'
        print self.numSpeziesVsGpcr

        pyplot.close()
        pyplot.xlabel('dataspace')
        pyplot.ylabel('rate')
        pyplot.xticks(range(4))
        pyplot.title('rate per dataspace')
        pyplot.bar(self.speziesVsGpcr, self.numSpeziesVsGpcr, color = '#BDBDBD')
        pyplot.savefig(os.path.join(stat6Dir, 'SpeziesVsGPCR.svg'), format='svg')

        # bestimmt den durchschnittlichen Score fuer jede Looplaenge
        self.scoreListY = []
        self.lenX = []

        gpcrFalse = 0
        gpcrTrue =0
        gpcrFoundList = {}

        for j, loop in enumerate(statSortedPdbLoopDictList) :
            if loop[9] == False :
                gpcrFalse +=1
            elif loop[9] == True :
                gpcrTrue +=1

        gpcrFoundList["LoopTotal"] = j+1
        gpcrFoundList["FalseTotal"] = gpcrFalse
        gpcrFoundList["TrueTotal"] = gpcrTrue

        with open (os.path.join(self.output_dir, "GPCR-Found.json"), 'w') as fp7 :
            json.dump(gpcrFoundList, fp7)


        for i in range(36) :
            # print i
            lengthList = [x for x in statSortedPdbLoopDictList if x[11] == i]
            lenList = len(lengthList)
            if lenList == 0 :
                continue
            # print lenList
            # print lengthList
            avgScore = sum([x[4] for x in lengthList], 0.0) / lenList
            self.scoreListY.append(avgScore)
            self.lenX.append(i)

        std = np.std(self.scoreListY)
        pyplot.xlabel('loop length')
        pyplot.ylabel('average score')
        pyplot.title('Average Score Per Loop Length')
        pyplot.bar(self.lenX, self.scoreListY, color = '#BDBDBD', yerr = std)
        pyplot.xticks(range(36))
        pyplot.savefig(os.path.join(stat1Dir, 'AvgScoreVsLooplengthPerLength.svg'),format='svg')

        self.scoreListY = []
        self.lenX = []

        for i in range(36) :
            # print i
            lengthList = [x for x in statSortedPdbLoopDictList if x[11] == i]
            lenList = len(lengthList)
            if lenList == 0 :
                continue
            # print lenList
            # print lengthList
            avgScore = sum([x[5] for x in lengthList], 0.0) / lenList
            self.scoreListY.append(avgScore)
            self.lenX.append(i)

        std = np.std(self.scoreListY)
        pyplot.close()
        pyplot.xlabel('loop length')
        pyplot.ylabel('average clashes')
        pyplot.title('Average Clashes Per Loop Length')
        pyplot.bar(self.lenX, self.scoreListY, color = '#BDBDBD', yerr = std)
        pyplot.xticks(range(36))
        pyplot.savefig(os.path.join(stat1Dir, 'AvgClashesVsLooplengthPerLength.svg'),format='svg')

        #bestimmt wie oft eine Erweiterung genutzt wurde, den Druchschnittsscore fuer jede Erweierung, Durchschnittliche Anzahl an Chlashes fuer jede Erweiterung
        self.rateY = []
        self.extensionX = []

        #print statSortedPdbLoopDictList

        for i in range(4) :
           extensionList = [x for x in statSortedPdbLoopDictList if x[0] == i]
           rateExtension = len(extensionList)
           self.rateY.append(rateExtension)
           self.extensionX.append(i)

        print 'Nummer of Extension:'
        print self.rateY

        std = np.std(self.rateY)
        pyplot.close()
        pyplot.xlabel('extension')
        pyplot.ylabel('rate')
        pyplot.xticks(range(4))
        pyplot.title('Rate per Extension')
        pyplot.bar(self.extensionX, self.rateY, color = '#BDBDBD', yerr = std)
        pyplot.savefig(os.path.join(stat2Dir, 'RateVsExtension.svg'), format='svg')

        self.rateY = []
        self.extensionX = []

        for i in range(4) :
            extensionList = [x for x in statSortedPdbLoopDictList if x[0] == i]
            if len(extensionList) != 0 :
                avgScoreExtension = sum([x[4] for x in extensionList]) / len(extensionList)
            else :
                avgScoreExtension = 0
            self.rateY.append(avgScoreExtension)
            self.extensionX.append(i)

        print 'Nummer of AvgScore:'
        print self.rateY

        std = np.std(self.rateY)
        pyplot.close()
        pyplot.xlabel('extension')
        pyplot.ylabel('avgScore')
        pyplot.xticks(range(4))
        pyplot.title('AvgScore per Extension')
        pyplot.bar(self.extensionX, self.rateY, color = '#BDBDBD', yerr = std)
        pyplot.savefig(os.path.join(stat2Dir, 'AvgScorePerExtension.svg'), format='svg')

        self.rateY = []
        self.extensionX = []

        for i in range(4) :
            extensionList = [x for x in statSortedPdbLoopDictList if x[0] == i]
            if len(extensionList) != 0 :
                avgScoreExtension = sum([x[5] for x in extensionList]) / len(extensionList)
            else :
                avgScoreExtension = 0
            self.rateY.append(avgScoreExtension)
            self.extensionX.append(i)

        print 'Nummer of AvgClashes:'
        print self.rateY

        std = np.std(self.rateY)
        pyplot.close()
        pyplot.xlabel('extension')
        pyplot.ylabel('avgClashes')
        pyplot.xticks(range(6))
        pyplot.title('AvgClashes per Extension')
        pyplot.bar(self.extensionX, self.rateY, color = '#BDBDBD', yerr = std)
        pyplot.savefig(os.path.join(stat2Dir, 'AvgClashesPerExtension.svg'), format='svg')

        #bestimmt fuer jede Erweiterung einzelt das Verhaeltnis von Score zuer Looplaenge
        self.extensionList0X = []
        self.extensionList0Y = []

        for i in range(4) :
            if i == 0 :
                extensionList0 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList0 = [x for x in extensionList0 if x[14] == k]
                    lenList0 = len(lengthList0)
                    if lenList0 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[4] for x in lengthList0], 0.0) / lenList0
                    self.extensionList0Y.append(avgScore)
                    self.extensionList0X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ScoreVsLoopLengthFor0,0-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std0 = np.std(self.extensionList0Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('score')
            pyplot.title('Score per Loop Length for 0,0')
            pyplot.bar(self.extensionList0X, self.extensionList0Y, color = '#BDBDBD', yerr = std0)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat3Dir, 'ScoreVsLoopLengthFor0,0.svg'),format='svg')


        self.extensionList1X = []
        self.extensionList1Y = []

        for i in range(4) :
            if i == 1 :
                extensionList1 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList1 = [x for x in extensionList1 if x[14] == k]
                    lenList1 = len(lengthList1)
                    if lenList1 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[4] for x in lengthList1], 0.0) / lenList1
                    self.extensionList1Y.append(avgScore)
                    self.extensionList1X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ScoreVsLoopLengthFor1,1-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std1 = np.std(self.extensionList1Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('score')
            pyplot.title('Score per Loop Length for 1,1')
            pyplot.bar(self.extensionList1X, self.extensionList1Y, color = '#BDBDBD', yerr = std1)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat3Dir, 'ScoreVsLoopLengthFor1,1.svg'), format='svg')

        self.extensionList2X = []
        self.extensionList2Y = []

        for i in range(4) :
            if i == 2 :
                extensionList2 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList2 = [x for x in extensionList2 if x[14] == k]
                    lenList2 = len(lengthList2)
                    if lenList2 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[4] for x in lengthList2], 0.0) / lenList2
                    self.extensionList2Y.append(avgScore)
                    self.extensionList2X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ScoreVsLoopLengthFor2,2-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std2 = np.std(self.extensionList2Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('score')
            pyplot.title('Score per Loop Length for 2,2')
            pyplot.bar(self.extensionList2X, self.extensionList2Y, color = '#BDBDBD', yerr = std2)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat3Dir, 'ScoreVsLoopLengthFor2,2.svg'),format='svg')

        self.extensionList3X = []
        self.extensionList3Y = []

        for i in range(4) :
            if i == 3 :
                extensionList3 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList3 = [x for x in extensionList3 if x[14] == k]
                    lenList3 = len(lengthList3)
                    if lenList3 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[4] for x in lengthList3], 0.0) / lenList3
                    self.extensionList3Y.append(avgScore)
                    self.extensionList3X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ScoreVsLoopLengthFor3,3-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std3 = np.std(self.extensionList3Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('score')
            pyplot.title('Score per Loop Length for 3,3')
            pyplot.bar(self.extensionList3X, self.extensionList3Y, color = '#BDBDBD', yerr = std3)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat3Dir, 'ScoreVsLoopLengthFor3,3.svg'),format='svg')

        barWidth = 0.2
        _, ax = pyplot.subplots()
        ex0 = ax.bar(np.array(self.extensionList0X) - barWidth, self.extensionList0Y, barWidth, color='#FA5882', yerr = std0, label = 'Extension_0')
        ex1 = ax.bar(np.array(self.extensionList1X) - 2*barWidth, self.extensionList1Y, barWidth, color='#5882FA', yerr = std1, label = 'Extension_1')
        ex2 = ax.bar(np.array(self.extensionList2X), self.extensionList2Y, barWidth, color='#FAAC58', yerr = std2, label = 'Extension_2')
        ex3 = ax.bar(np.array(self.extensionList3X) + barWidth, self.extensionList3Y, barWidth, color='#58FA58', yerr = std3, label = 'Extension_3')
        pyplot.xticks(range(36))
        ax.set_xlabel('ori loop length')
        ax.set_ylabel('score')
        ax.set_title('Score per Loop Length')
        pyplot.legend()
        #ax.set_xticks([x + barWidth for x in range(41)])
        pyplot.tight_layout()
        pyplot.savefig(os.path.join(stat3Dir, 'ScoreVsLoopLength.svg'),format='svg')
        #bestimmt durschnittlichen Score des durchschnittlichen Scores von und fuer pdb 1-10

        self.extensionList4X = []
        self.extensionList4Y = []

        # for i in range(6) :
        #     if i == 4 :
        #         extensionList4 = [x for x in statSortedPdbLoopDictList if x[0] == i]
        #         for k in range(36) :
        #             # print i
        #             lengthList4 = [x for x in extensionList4 if x[14] == k]
        #             lenList4 = len(lengthList4)
        #             if lenList4 == 0 :
        #                 continue
        #             # print lenList
        #             # print lengthList
        #             avgScore = sum([x[4] for x in lengthList4], 0.0) / lenList4
        #             self.extensionList4Y.append(avgScore)
        #             self.extensionList4X.append(k)
        #
        # if len(self.extensionList0Y) == 0 :
        #     with open (os.path.join(stat3Dir, "ScoreVsLoopLengthFor4,4-notExisting.txt"),'w' ) as fp10 :
        #         fp10.write('wurde nicht ausgewaehlt')
        # else :
        #     std4 = np.std(self.extensionList4Y)
        #     pyplot.close()
        #     pyplot.xlabel('ori loop length')
        #     pyplot.ylabel('score')
        #     pyplot.title('Score per Loop Length for 4,4')
        #     pyplot.bar(self.extensionList4X, self.extensionList4Y, color = '#BDBDBD', yerr = std4)
        #     pyplot.xticks(range(36))
        #     pyplot.savefig(os.path.join(stat3Dir, 'ScoreVsLoopLengthFor4,4.svg'),format='svg')
        #
        # self.extensionList5X = []
        # self.extensionList5Y = []

        # for i in range(6) :
        #     if i == 5 :
        #         extensionList5 = [x for x in statSortedPdbLoopDictList if x[0] == i]
        #         for k in range(36) :
        #             # print i
        #             lengthList5 = [x for x in extensionList5 if x[14] == k]
        #             lenList5 = len(lengthList5)
        #             if lenList5 == 0 :
        #                 continue
        #             # print lenList
        #             # print lengthList
        #             avgScore = sum([x[4] for x in lengthList5], 0.0) / lenList5
        #             self.extensionList5Y.append(avgScore)
        #             self.extensionList5X.append(k)
        #
        # if len(self.extensionList0Y) == 0 :
        #     with open (os.path.join(stat3Dir, "ScoreVsLoopLengthFor5,5-notExisting.txt"),'w' ) as fp10 :
        #         fp10.write('wurde nicht ausgewaehlt')
        # else :
        #     std5 = np.std(self.extensionList5Y)
        #     pyplot.close()
        #     pyplot.xlabel('ori loop length')
        #     pyplot.ylabel('score')
        #     pyplot.title('Score per Loop Length for 5,5')
        #     pyplot.bar(self.extensionList5X, self.extensionList5Y, color = '#BDBDBD', yerr = std5)
        #     pyplot.xticks(range(36))
        #     pyplot.savefig(os.path.join(stat3Dir, 'ScoreVsLoopLengthFor5,5.svg'),format='svg')

        #bestimmt fuer jede Erweiterung einzelt das Verhaeltnis von Clashes zur Looplaenge
        self.extensionList0X = []
        self.extensionList0Y = []

        for i in range(4) :
            if i == 0 :
                extensionList0 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList0 = [x for x in extensionList0 if x[14] == k]
                    lenList0 = len(lengthList0)
                    if lenList0 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[5] for x in lengthList0], 0.0) / lenList0
                    self.extensionList0Y.append(avgScore)
                    self.extensionList0X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ClashVsLoopLengthFor0,0-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std0 = np.std(self.extensionList0Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('clash')
            pyplot.title('Clashes per Loop Length for 0,0')
            pyplot.bar(self.extensionList0X, self.extensionList0Y, color = '#BDBDBD', yerr = std0)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat4Dir, 'ClashVsLoopLengthFor0,0.svg'), format='svg')

        self.extensionList1X = []
        self.extensionList1Y = []

        for i in range(4) :
            if i == 1 :
                extensionList1 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList1 = [x for x in extensionList1 if x[14] == k]
                    lenList1 = len(lengthList1)
                    if lenList1 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[5] for x in lengthList1], 0.0) / lenList1
                    self.extensionList1Y.append(avgScore)
                    self.extensionList1X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ClashVsLoopLengthFor1,1-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std1 = np.std(self.extensionList1Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('clash')
            pyplot.title('Clashes per Loop Length for 1,1')
            pyplot.bar(self.extensionList1X, self.extensionList1Y, color = '#BDBDBD', yerr = std1)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat4Dir, 'ClashVsLoopLengthFor1,1.svg'),format='svg')

        self.extensionList2X = []
        self.extensionList2Y = []

        for i in range(4) :
            if i == 2 :
                extensionList2 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList2 = [x for x in extensionList2 if x[14] == k]
                    lenList2 = len(lengthList2)
                    if lenList2 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[5] for x in lengthList2], 0.0) / lenList2
                    self.extensionList2Y.append(avgScore)
                    self.extensionList2X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ClashVsLoopLengthFor2,2-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std2 = np.std(self.extensionList2Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('clash')
            pyplot.title('Clashes per Loop Length for 2,2')
            pyplot.bar(self.extensionList2X, self.extensionList2Y, color = '#BDBDBD', yerr = std2)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat4Dir, 'ClashVsLoopLengthFor2,2.svg'),format='svg')

        self.extensionList3X = []
        self.extensionList3Y = []

        for i in range(4) :
            if i == 3 :
                extensionList3 = [x for x in statSortedPdbLoopDictList if x[0] == i]
                for k in range(36) :
                    # print i
                    lengthList3 = [x for x in extensionList3 if x[14] == k]
                    lenList3 = len(lengthList3)
                    if lenList3 == 0 :
                        continue
                    # print lenList
                    # print lengthList
                    avgScore = sum([x[5] for x in lengthList3], 0.0) / lenList3
                    self.extensionList3Y.append(avgScore)
                    self.extensionList3X.append(k)

        if len(self.extensionList0Y) == 0 :
            with open (os.path.join(stat3Dir, "ClashVsLoopLengthFor3,3-notExisting.txt"),'w' ) as fp10 :
                fp10.write('wurde nicht ausgewaehlt')
        else :
            std3 = np.std(self.extensionList3Y)
            pyplot.close()
            pyplot.xlabel('ori loop length')
            pyplot.ylabel('clash')
            pyplot.title('Clashes per Loop Length for 3,3')
            pyplot.bar(self.extensionList3X, self.extensionList3Y, color = '#BDBDBD', yerr = std3)
            pyplot.xticks(range(36))
            pyplot.savefig(os.path.join(stat4Dir, 'ClashVsLoopLengthFor3,3.svg'),format='svg')

        pyplot.close()

        barWidth = 0.2
        _, ax = pyplot.subplots()
        ex0 = ax.bar(np.array(self.extensionList0X) - 2*barWidth, self.extensionList0Y, barWidth, color='#FA5882', yerr = std0, label = 'Extension_0')
        ex1 = ax.bar(np.array(self.extensionList1X) - barWidth, self.extensionList1Y, barWidth, color='#5882FA', yerr = std1, label = 'Extension_1')
        ex2 = ax.bar(np.array(self.extensionList2X), self.extensionList2Y, barWidth, color='#FAAC58', yerr = std2, label = 'Extension_2')
        ex3 = ax.bar(np.array(self.extensionList3X) + barWidth, self.extensionList3Y, barWidth, color='#58FA58', yerr = std3, label = 'Extension_3')
        pyplot.xticks(range(36))
        ax.set_xlabel('loop length')
        ax.set_ylabel('clash')
        ax.set_title('Clash per Loop Length')
        pyplot.legend()
        #ax.set_xticks([x + barWidth for x in range(41)])
        pyplot.tight_layout()
        pyplot.savefig(os.path.join(stat4Dir, 'ClashVsLoopLength.svg'),format='svg')


        # self.extensionList4X = []
        # self.extensionList4Y = []
        #
        # for i in range(6) :
        #     if i == 4 :
        #         extensionList4 = [x for x in statSortedPdbLoopDictList if x[0] == i]
        #         for k in range(36) :
        #             # print i
        #             lengthList4 = [x for x in extensionList4 if x[14] == k]
        #             lenList4 = len(lengthList4)
        #             if lenList4 == 0 :
        #                 continue
        #             # print lenList
        #             # print lengthList
        #             avgScore = sum([x[5] for x in lengthList4], 0.0) / lenList4
        #             self.extensionList4Y.append(avgScore)
        #             self.extensionList4X.append(k)
        #
        # if len(self.extensionList0Y) == 0 :
        #     with open (os.path.join(stat3Dir, "ClashVsLoopLengthFor4,4-notExisting.txt"),'w' ) as fp10 :
        #         fp10.write('wurde nicht ausgewaehlt')
        # else :
        #     std4 = np.std(self.extensionList4Y)
        #     pyplot.close()
        #     pyplot.xlabel('ori loop length')
        #     pyplot.ylabel('clash')
        #     pyplot.title('Clashes per Loop Length for 4,4')
        #     pyplot.bar(self.extensionList4X, self.extensionList4Y, color = '#BDBDBD', yerr = std4)
        #     pyplot.xticks(range(36))
        #     pyplot.savefig(os.path.join(stat4Dir, 'ClashVsLoopLengthFor4,4.svg'), format='svg')
        #
        # self.extensionList5X = []
        # self.extensionList5Y = []
        #
        # for i in range(6) :
        #     if i == 5 :
        #         extensionList5 = [x for x in statSortedPdbLoopDictList if x[0] == i]
        #         for k in range(36) :
        #             # print i
        #             lengthList5 = [x for x in extensionList5 if x[14] == k]
        #             lenList5 = len(lengthList5)
        #             if lenList5 == 0 :
        #                 continue
        #             # print lenList
        #             # print lengthList
        #             avgScore = sum([x[5] for x in lengthList5], 0.0) / lenList5
        #             self.extensionList5Y.append(avgScore)
        #             self.extensionList5X.append(k)
        #
        # if len(self.extensionList0Y) == 0 :
        #     with open (os.path.join(stat3Dir, "ClashVsLoopLengthFor5,5-notExisting.txt"),'w' ) as fp10 :
        #         fp10.write('wurde nicht ausgewaehlt')
        # else :
        #     std5 = np.std(self.extensionList5Y)
        #     pyplot.close()
        #     pyplot.xlabel('ori loop length')
        #     pyplot.ylabel('clash')
        #     pyplot.title('Clashes per Loop Length for 5,5')
        #     pyplot.bar(self.extensionList5X, self.extensionList5Y, color = '#BDBDBD', yerr = std5)
        #     pyplot.xticks(range(36))
        #     pyplot.savefig(os.path.join(stat4Dir, 'ClashVsLoopLengthFor5,5.svg'), format='svg')

        # self.avgScoreY = []
        # self.modelNumX = []
        # for i in range(11) :
        #     avgList = [x for x in statAvgScoreList if x [0] == i]
        #     lenAvg = len(avgList)
        #     if lenAvg == 0 :
        #         continue
        #     avgScore = sum(x[1] for x in avgList)/ lenAvg
        #     self.avgScoreY.append(avgScore)
        #     std = np.std(self.avgScoreY)
        #     self.modelNumX.append(i)
        #
        #     # print self.avgScoreY
        #     # print '==========='
        #     # print self.modelNumX
        #
        # pyplot.close()
        # pyplot.xlabel('model')
        # pyplot.ylabel('avgScore')
        # pyplot.title('AvgScore per Model')
        # pyplot.bar(self.modelNumX, self.avgScoreY, color = '#BDBDBD', yerr = std)
        # pyplot.savefig(os.path.join(stat5Dir, 'AvgScore per Model'))

# fuer Failed.txt zusammenfassen
class SSFEFailed (PyTool, ProviMixin ):
    args = [
        _("out_dataSet", type="dir")
    ]

    def func ( self ) :

        zipf = zipfile.ZipFile( os.path.join(self.output_dir, "failed.zip"), 'w', zipfile.ZIP_DEFLATED)
        include = set(['Failed.txt', 'Failed1.txt', 'Failed2.txt', 'Failed3.txt'])
        exclude = ['link_it_0,0', 'link_it_1,1', 'link_it_2,2', 'link_it_3,3', 'Helix8']
        for root, dirs, files in os.walk(self.out_dataSet, topdown=True) :
            dirs[:] = [d for d in dirs if d not in exclude]
            if not set(files).intersection(include) and 'Result_Table.csv' in files:
                continue
            for file in files :
                zipf.write(os.path.join(os.path.relpath(root), file))
        zipf.close()



# fuert SSFE fuer einen Datensatz aus
class SSFEZip (PyTool, ProviMixin ):
    args = [
        _("out_dataSet", type="dir")
    ]

    def func ( self ) :

        # path = os.path.join(outJobDir)

        zipf = zipfile.ZipFile( os.path.join(self.output_dir, "results.zip"), 'w', zipfile.ZIP_DEFLATED)
        exclude = ['link_it_0,0', 'link_it_1,1', 'link_it_2,2', 'link_it_3,3', 'Helix8']
        # 'link_it_4,4','link_it_5,5',
        for root, dirs, files in os.walk(self.out_dataSet, topdown=True) :
            dirs[:] = [d for d in dirs if d not in exclude]
            for file in files :
                zipf.write(os.path.join(os.path.relpath(root), file))
        zipf.close()


class SSFEMultiLinkIt( PyTool, ProviMixin ):
    args = [
        _("loop_jobs", type="dir"),
        _( "GPCRscore", type="int", default=1 ),
        _( "Speciesscore", type="int", default=2 ),
        _( "clashOut", type="float", default=0.75 ),
        _( "numclashes", type="int", default=10)
    ]

    def func ( self ) :

        os.chdir(self.output_dir)

        self.loop_jobs = os.path.abspath(self.loop_jobs)

        #enthaelt Liste mit allen Dateien in denen die Proteinen + ihr Loops stehen
        jobDirs = os.listdir (self.loop_jobs)

        for jobDir in jobDirs :
            jobDirPath = os.path.join(self.loop_jobs, jobDir)
            jobFileList = {}
            pdbFileList = []

            gzipFiles = [f for f in os.listdir(jobDirPath) if f.endswith('.gz')]

            os.chdir(jobDirPath)
            #entpacke eingabe daten
            for gzipFile in gzipFiles :
                with gzip.open(gzipFile, 'rb') as f:
                    with open(os.path.splitext(gzipFile)[0], 'w') as outfile:
                        for line in f:
                            outfile.write(line)
            os.chdir('..')

            inputFiles = os.listdir(jobDirPath)

            for inputFile in inputFiles :
                if  (os.path.splitext(inputFile)[1] == '.pdb') :
                    pdbFileList.append(inputFile)
                else :
                    nrjob = re.compile("([a-zA-Z_]+)([0-9]+)")
                    #print nrjob
                    numberMatch = nrjob.match(inputFile.split("_")[-1])
                    if numberMatch :
                        numberjob = numberMatch.group(2)
                        jobFileList[numberjob]=inputFile

            # print pdbFileList
            # print jobFile

            outJobDir = self.subdir(jobDir + '_Result')

            os.chdir(outJobDir)
            #print pdbFileList
            #print '==============================================='
            for pdbFile in pdbFileList :
                outJobFileDir = self.subdir(os.path.join(outJobDir, os.path.splitext(pdbFile)[0]))

                os.chdir(outJobFileDir)

                # print '|',pdbFileList, '|',jobFile,'|'
                #print jobDirPath, pdbFile
                nr = re.compile("([a-zA-Z_]+)([0-9]+)")
                number = nr.match(pdbFile).group(2)
                jobFile = jobFileList[number]
                #print jobFile, pdbFile
                #print '==============================================='
                SSFELinkIt(os.path.join(jobDirPath, pdbFile), os.path.join(jobDirPath, jobFile), Speciesscore=self.Speciesscore, GPCRscore=self.GPCRscore, verbose=self.verbose, debug=True)

                os.chdir('..')

            os.chdir('..')

        # path = os.path.join(outJobDir)
        #
        #
        # zipf = zipfile.ZipFile( os.path.join(self.output_dir, "results.zip"), 'w', zipfile.ZIP_DEFLATED)
        # exclude = ['link_it_0,0', 'link_it_1,1', 'link_it_2,2', 'link_it_3,3']
        # for root, dirs, files in os.walk(path, topdown=True) :
        #     dirs[:] = [d for d in dirs if d not in exclude]
        #     for file in files :
        #         zipf.write(os.path.join(root, file))
        # zipf.close()



class SSFELinkIt( PyTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "loop_file", type="file", ext="txt" ),
        _( "extension", type="list", nargs=2, action="append",
           help="int,int", default=[3,3] ),
        _( "GPCRscore", type="int", default=5 ),
        _( "Speciesscore", type="int", default=5 ),
        _( "clashOut", type="float", default=0.75 ),
        _( "numclashes", type="int", default=10)

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
                    # print line[11:17].split()
                    if line[11:17].split()[0] != "OXT":
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
        #fur insertLoop --> einbau aller loops in ein pdb.file
        self.loopPosList = []
        self.loopBracketAminosHelix8 = {}


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
                    seqLength = 35
                    if len(seq) > seqLength :
                        # self.loopTasks = []
                        self.loopPosList.append([0, 0])
                        with open(os.path.join( self.output_dir, "ERROR_LoopSequenzToLong.txt"), 'w') :
                            pass
                        continue
                        # loopBrackets = []
                        # loopBracketAminos = {}
                        # with open(os.path.join( self.output_dir, "ERROR_LoopSequenzToLong.txt"), 'w') :
                        #     pass
                        # return
                    # fehlerausgabe fuer falsche eingaben
                    #or (int(start) < int(end) and len(seq) != 0 and not line[0:7] == 'between' )
                    if int(start) > int(end) and len(seq) != 0 or (int(start) < int(end) and len(seq) == 0 and not line[0:7] == 'between' ) :
                        # print int(start)
                        # print int(end)
                        # print len(seq)
                        # print '====================='
                        self.loopTasks = []
                        loopBrackets = []
                        loopBracketAminos = {}
                        with open(os.path.join( self.output_dir, "ERROR_InvalidInput.txt"), 'w') :
                            pass
                        return

                    # loop muss berechnet werden, damit in index keine lueck, wird unabhaenig der andern loops gemacht und nur mit Erweiterung 3,3
                    # if start > end and len(seq) == 0 or :
                    #     startInt = int(end)
                    #     endInt = int(start)
                    # 'TM7_H8'

                    if line[0:7] == 'between' and len(seq) == 0 :
                        startInt = int(start)
                        #print start
                        endInt = int(end)
                        #print end
                        self.loopBracketAminosHelix8 = {}#dict.fromkeys([startInt-1, startInt, endInt,endInt+1])
                        #print self.loopBracketAminosHelix8
                        self.loopPosList.append([0, 0])
                        continue

                    endInt2 = int(end)+1
                    with open( self.pdb_file, 'r' ) as fp :
                        for line in fp :
                            if line.startswith("ATOM") :
                                position = int(line[22:26])
                                last_res=position
                    # print last_res
                    # print '------------------'
                    # print endInt2+3
                    if last_res < endInt2+3 :
                        with open ('ERROR_TMH7_H8_2.txt' , 'w') as fp7 :
                            fp7.write( 'TMH7/H8 zu nah am Ende des Rezeptors - keine Berechnung moeglich')
                        continue

                    self.loopCount += 1

                    #Eingabeliste aller zu bearbeitenden Loops aus Eingabedatei
                    self.loopTasks.append([int(start)-1 , int(end)+1 , seq])
                    self.loopPosList.append([start, end])

                    #List von Start- und Endwertern der erweiterten Loops
                    startInt = int(start)-1
                    endInt = int(end)+1
                    # with open( self.pdb_file, 'r' ) as fp :
                    #     for line in fp :
                    #         if line.startswith("ATOM") :
                    #             position = int(line[22:26])
                    #             last_res=position
                    # print last_res
                    # print '------------------'
                    # print endInt+3
                    # if last_res < endInt+3 :
                    #     continue
                    loopBrackets.append(([startInt-2,startInt-1,startInt], [endInt, endInt+1, endInt+2]))
                    # startInt-4,startInt-3, , endInt+3, endInt+4
                    loopBracketAminos.update(dict.fromkeys(loopBrackets[-1][0]))
                    loopBracketAminos.update(dict.fromkeys(loopBrackets[-1][1]))

                    # for key in loopBracketAminos :
                    #     for loop in loopBracketAminos[key] :
                    #         for pos in loopBracketAminos[key][loop] :
                    #             if pos > last_res :
                    #
                    #             else :




        self.loopPosList += (7 * [[0, 0]])
        #print self.loopBracketAminosHelix8


        #print loopBracketAminos
        #bestimmt Listen mit N und C Terminus der einzelnen Loops
        last_res = 9999
        with open( self.pdb_file, 'r' ) as fp :
            for line in fp :
                if line.startswith("ATOM") :
                    position = int(line[22:26])
                    last_res=position
                    if position in loopBracketAminos :
                        loopBracketAminos[position] = line[17:20]

        if None in loopBracketAminos.values() :
            with open ('Failed.txt' , 'w') as fp :
                fp.write( 'Loopsequenze und Template stimmen nicht ueberein. Looppositionen im Template ueberpruefen')
            return

        #ersetzt 3-Buchstabencode durch 1-Buchstabencode
        for position, aminoTriplet in loopBracketAminos.iteritems() :
            # print position,aminoTriplet
            loopBracketAminos[position] = SSFELinkIt.AminoDict[aminoTriplet]


        #fuer loop helix8
        if not self.loopBracketAminosHelix8 == {} :
            failed = False
            with open( self.pdb_file, 'r' ) as fp :
                for line in fp :
                    if line.startswith("ATOM") :
                        position = int(line[22:26])
                        last_res=position
                        if position in self.loopBracketAminosHelix8 :
                            self.loopBracketAminosHelix8[position] = line[17:20]

            for position, aminoTriplet in self.loopBracketAminosHelix8.iteritems() :
                #print position,aminoTriplet
                #print self.loopBracketAminosHelix8
                if last_res >= max(self.loopBracketAminosHelix8.keys()):
                    self.loopBracketAminosHelix8[position] = SSFELinkIt.AminoDict[aminoTriplet]
                else :
                    with open ('ERROR_TMH7_H8_1.txt' , 'a') as fp :
                        fp.write( 'Luecke zwischen TMH7 und Helix8 konnte fuer die Visualisierung nicht geschlossen werden. Sequenzpositionen nicht korrekt/verarbeitbar.')
                    failed = True
            if not failed :


                #print self.loopBracketAminosHelix8

                #print loopBracketAminos
                #-1 und +1 fuer Gly/ wir spaeter wieder abgeschnitten
                startLoop = min(self.loopBracketAminosHelix8.keys()) -1
                #print type(last_res)
                seqLoop = ''
                for key in sorted(self.loopBracketAminosHelix8) :
                    seqLoop += self.loopBracketAminosHelix8[key]
                if max(self.loopBracketAminosHelix8.keys()) == last_res:
                    endLoop = max(self.loopBracketAminosHelix8.keys())
                    seqLoop = seqLoop[:-1]
                else:
                    endLoop = max(self.loopBracketAminosHelix8.keys()) +1


                # c = chain
                loopTasksListHelix8 = [[str(startLoop) + ':' + c, str(endLoop) + ':' + c, seqLoop]]
                #print self.loopBracketAminosHelix8
                # print c
                # print loopTasksListHelix8
                # print self.loopBracketAminosHelix8

                rrrrr = MultiLinkIt(self.ori_file, input=loopTasksListHelix8,
                        **copy_dict( kwargs, run=False, gpcrDB=False, 
                                    output_dir=self.subdir("Helix8" ), verbose=True, debug=True ))
                rrrrr()
                
                loopFile = []
                loopLineList = []
                with open(self.subdir("Helix8/link_it_0") + "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker2.pdb", 'r') as fp :
                    loopFile = fp.readlines()

                modelBegin = -1
                modelEnd = -1
                breakout = False
                for i, line in enumerate(loopFile) :
                    if i==0 and line.startswith("END"):
                        with open ('FailedTH7_H8_2.txt' , 'a') as fp :
                            fp.write( 'Luecke zwischen TH7 und Helix8 konnte fuer die Visualisierung nicht geschlossen werden.')
                        breakout = True
                    else:
                        if line[0:5] == 'MODEL' and line.split()[1] == '1' :
                            modelBegin = i
                        if line[0:6] == 'ENDMDL' and modelBegin != -1 :
                            modelEnd = i
                            break
                if not breakout:

                    loopLineList.append(loopFile[modelBegin + 5 : modelEnd - 4])

                    # print modelBegin
                    # print modelEnd
                    # print loopLineList


                    new_pdb_file = (os.path.splitext(self.ori_file)[0] + 'Helix8.pdb')
                    with open( self.ori_file, 'r') as fp, open (new_pdb_file, 'w') as fp2 :
                        for line in fp :
                            #Filter ANISOU-Zeilen raus, falls vorhanden und kopiert alle anderen Zeilen
                            if line[0:6] != 'ANISOU' :
                                fp2.write(line)
                    for loop in loopLineList :
                        #print loop
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



                        shutil.move(new_pdb_file + '.tmp', self.ori_file)
                        
        #
        # fuer restliche loops
        self.oriSeqDict = {key: [] for key in ([0,1,2,3])}
        # ,4,5
        for i in range (4) :
            loopTasksList = []
            for start,end,seq in self.loopTasks :
                for j in range(i) :
                    seq = loopBracketAminos[start - j ] + seq + loopBracketAminos[end + j]
                loopTasksList.append([str(start - i) + ':' + c, str(end + i) + ':' + c, seq])
            # print '------------'
            # print loopTasksList
            for k in loopTasksList :
                self.oriSeqDict[i].append(k[2])



            self.multiLinkIts.append(MultiLinkIt(self.ori_file, input = loopTasksList,
                    **copy_dict( kwargs, run=False,
                                output_dir=self.subdir("membranDB/link_it_%i,%i" % (i,i) ), gpcrDB=False, last_res = last_res, verbose=True, debug=True )))

            self.multiLinkIts.append(MultiLinkIt(self.ori_file, input = loopTasksList,
                    **copy_dict( kwargs, run=False,
                                output_dir=self.subdir("GPCRDB/link_it_%i,%i" % (i,i) ), gpcrDB=True, last_res = last_res, verbose=True, debug=True )))

        #print self.oriSeqDict

    def compareLoops( self, x,y):
            if x[0] > y[0] :
                if (x[0]*self.clashOut) < y[0] :
                    if x[1] <= y[1] :
                        return 1
                    else :
                        return -1
                else :
                    return 1
            elif x[0] < y[0] :
                if x[0] > (y[0]*self.clashOut) :
                    if x[1] < y[1] :
                        return 1
                    else :
                        return -1
                else :
                    return -1
            else :
                if x[1] < y[1] :
                    return 1
                else :
                    return -1
                return 0

    def func( self ):
        if len(self.multiLinkIts) == 0 :
            return

        for m in self.multiLinkIts :
            m()


        # print self.loopCount

        self.gpcrListe = ["1EDN", "1F88", "1GZM", "1HZX", "1JFP", "1L9H", "1LN6", "1U19", "1V6R", "2G87",
                          "2HPY", "2I35", "2I36", "2I37", "2J4Y", "2LNL", "2LOT", "2LOU", "2LOV", "2LOW",
                          "2PED", "2R4R", "2R4S", "2RH1", "2VT4", "2X72", "2Y00", "2Y01", "2Y02", "2Y03",
                          "2Y04", "2YCW", "2YCX", "2YCY", "2YCZ", "2YDO", "2YDV", "2Z73", "2ZIY", "3AYM",
                          "3AYN", "3C9L", "3C9M", "3CAP", "3D4S", "3DQB", "3EML", "3KJ6", "3NY8", "3NY9",
                          "3NYA", "3OAX", "3ODU", "3OE0", "3OE6", "3OE8", "3OE9", "3P0G", "3PBL", "3PDS",
                          "3PQR", "3PWH", "3PXO", "3QAK", "3REY", "3RFM", "3RZE", "3SN6", "3UON", "3UZA",
                          "3UZC", "3V2W", "3V2Y", "3VG9", "3VGA", "3VW7", "3ZEV", "3ZPQ", "3ZPR", "4A4M",
                          "4AMI", "4AMJ", "4BEY", "4BEZ", "4BUO", "4BV0", "4BVN", "4BWB", "4DAJ", "4DJH",
                          "4DKL", "4EA3", "4EIY", "4EJ4", "4GBR", "4GPO", "4GRV", "4IAQ", "4IAR", "4IB4",
                          "4J4Q", "4JKV", "4K5Y", "4L6R", "4MBS", "4MQE", "4MQF", "4MQS", "4MQT", "4MR7",
                          "4MR8", "4MR9", "4MRM", "4MS1", "4MS3", "4MS4", "4N4W", "4N6H", "4NC3", "4NTJ",
                          "4OO9", "4OR2", "4PHU", "4PXF", "4PXZ", "4PY0", "4RWA", "4RWD", "4RWS", "4S0V",
                          "4U14", "4UG2", "4UHR", "4WW3", "4X1H", "4XEE", "4XES", "4XNV", "4XNW", "4XT1",
                          "4XT3", "4YAY", "4Z34", "4Z35", "4Z36", "4Z9G", "4ZJ8", "4ZJC", "4ZUD", "4ZWJ",
                          "5A8E", "5C1M", "5CGC", "5CGD", "5CXV", "5D5A", "5D5B", "5D6L", "5DGY", "5DHG",
                          "5DHH", "5DSG", "5DYS", "5EN0", "5F8U", "5G53", "5GLH", "5IU4", "5IU7", "5IU8",
                          "5IUA", "5IUB", "5K2A", "5K2B", "5K2C", "5K2D", "5L7D", "5L7I", "5LWE", "5T04",
                          "5T1A", "5TGZ", "5TVN", "5U09"]

        self.gpcrDict = {'X_5HT2B': ["4IB4", "4NC3", "5TVN"], 'X_ENLYS': ["2RH1"],
                         'X_GABR1': ["4MQE", "4MQF", "4MR7", "4MR8", "4MR9", "4MRM", "4MS1", "4MS3", "4MS4"],
                         'X_NTR1': ["3ZEV", "4BUO", "4BV0", "4BWB", "4GRV", "4XEE", "4XES", "5T04"],
                         'X_GABR2': ["4MQE", "4MQF", "4MR8", "4MR9", "4MRM", "4MS1", "4MS3", "4MS4"],
                         'X_SMO': ["4JKV", "4N4W", "5L7D", "5L7I"], 'X_OX2R': ["4S0V"], 'X_CCR2': ["5T1A"],
                         'X_CCR5': ["4MBS"], 'X_CCR9': ["5LWE"], 'X_GLR': ["4L6R"], 'X_EDNRB': ["5GLH"],
                         'X_OPRK': ["4DJH"], 'X_P2Y12': ["4NTJ", "4PXZ", "4PY0"], 'X_P2RY1': ["4XNV", "4XNW"],
                         'X_OPRM': ["4DKL", "5C1M"], 'X_CRFR1': ["4K5Y", "4Z9G"], 'X_S1PR1': ["3V2W", "3V2Y"],
                         'X_OPRD': ["4EJ4", "4N6H", "4RWA", "4RWD"], 'X_CXCR1': ["2LNL"], 'X_OPRX': ["4EA3", "5DHG", "5DHH"],
                         'X_HRH1': ["3RZE"], 'X_EDN1': ["1EDN", "1V6R"], 'X_GRM5': ["4OO9", "5CGC", "5CGD"], 'X_GRM1': ["4OR2"],
                         'X_CNR1': ["5TGZ", "5U09"],'X_AA2AR': ["2YDO", "2YDV", "3EML", "3PWH", "3QAK", "3REY", "3RFM", "3UZA", "3UZC", "3VG9", "3VGA", "4EIY",
                                                            "4UG2", "4UHR", "5G53", "5IU4", "5IU7", "5IU8", "5IUA", "5IUB", "5K2A", "5K2B", "5K2C", "5K2D"],
                         'X_5HT1B': ["4IAQ", "4IAR"], 'X_CXCR4': ["3ODU", "3OE0", "3OE6", "3OE8", "3OE9", "4RWS"], 'X_ACM3': ["4DAJ", "4U14"],
                         'X_ACM2': ["3UON", "4MQS", "4MQT"], 'X_ACM1': ["5CXV"], 'X_ACM4': ["5DSG"], 'X_FFAR1': ["4PHU"], 'X_OX1R': ["4ZJ8", "4ZJC"],
                         'X_US28': ["4XT1", "4XT3"], 'X_OPSD': ["1F88", "1GZM", "1HZX", "1JFP", "1L9H", "1LN6", "1U19", "2G87", "2HPY", "2I35", "2I36", "2I37",
                                                            "2J4Y", "2PED", "2X72", "2Z73", "2ZIY", "3AYM", "3AYN", "3C9L", "3C9M", "3CAP", "3DQB", "3OAX",
                                                            "3PQR", "3PXO", "4A4M", "4BEY", "4BEZ", "4J4Q", "4PXF", "4WW3", "4X1H", "4ZWJ", "5DGY", "5DYS", "5EN0"],
                         'X_DRD3': ["3PBL"], 'X_APJ': ["2LOT", "2LOU", "2LOV", "2LOW"], 'X_ADRB1': ["2VT4", "2Y00", "2Y01", "2Y02", "2Y03", "2Y04", "2YCW",
                                                                                              "2YCX", "2YCY", "2YCZ", "3ZPQ", "3ZPR", "4AMI", "4AMJ",
                                                                                              "4BVN", "4GPO", "5A8E", "5F8U"],
                         'X_ADRB2': ["2R4R", "2R4S", "3D4S", "3KJ6", "3NY8", "3NY9", "3NYA", "3P0G", "3PDS", "3SN6", "4GBR", "5D5A", "5D5B", "5D6L"],
                         'X_PAR1': ["3VW7"], 'X_AGTR1': ["4YAY", "4ZUD"], 'X_LPAR1': ["4Z34", "4Z35", "4Z36"]}



        #erstellt Dict mit den 6 besten berechneten Loops pro Loop
        #list, damit keine tamplate doppelt vorkommt (z.b.4z36)

        speciesName = os.path.splitext(os.path.basename(self.loop_file))[0].split("_")[0]
        #print "speciesName", speciesName

        def firstsort( outputDir, memDB=False ):
            #loopDict enthaelt fuer jeden Loop die 3 besten berechneten Loops pro Erweiterung:(0,0);(1,1);(2,2);(3,3)
            temp_loopDict = {}
            #zum Laden der einzelnen Loop-json zur Auswahl der 3 besten Loops nach Score und Anzahl Clashes
            singleLoopDict = {}
            templateList = []
            for i in range(4) :
                for j in range(self.loopCount) :
                    jsonPfad = os.path.join( outputDir,"link_it_%i,%i/link_it_%i" % (i,i,j) )+ "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker.json"
                    if os.path.isfile(jsonPfad) :
                        with open(os.path.join( outputDir,"link_it_%i,%i/link_it_%i" % (i,i,j) )+ "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker.json", 'r') as fp :
                            singleLoopDict = json.load(fp)
                            loop = temp_loopDict.get(j, [])
                            #print type(loopDict)
                            #print (os.path.join( self.output_dir,"link_it_%i,%i/link_it_%i" % (i,i,j) )+ "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker.json")

                        sortedSingleLoopDict = []
                            # damit man noch weiss welches Model(1-100) zu den Werten gehoert
                        for key in singleLoopDict["linker"] :
                            sortedSingleLoopDict += [[int(key)] + singleLoopDict["linker"][key]]
                            # print k
                            # print '========'

                        # print '========================================'
                        #print sortedSingleLoopDict
                        #loecht eintraege wenn tamplate doppelt vorkommt (z.b.4z36)
                        # for loopElement in sortedSingleLoopDict :
                        #     if loopElement[4] in templateList :
                        #         sortedSingleLoopDict.remove(loopElement)
                        #     else :
                        #         templateList.append(loopElement[4])

                        # Sezies hoeher werten
                        for index, loopEntrie in enumerate(sortedSingleLoopDict) :
                            try:
                                if loopEntrie[4].upper() in self.gpcrDict["X_" + speciesName] :
                                    # if loopEntrie[4].upper() == '3ODU' :
                                    #     print loopEntrie
                                    sortedSingleLoopDict[index][2] *= self.Speciesscore
                                    sortedSingleLoopDict[index].append(True)
                                else:
                                    sortedSingleLoopDict[index].append(False)
                            except:
                                sortedSingleLoopDict[index].append(False)

                        # print sortedSingleLoopDict
                        # break
                        # GPCR Score bei gefunden verdoppelt
                        for index, loopEntrie in enumerate(sortedSingleLoopDict) :
                            if loopEntrie[4].upper() in self.gpcrListe :
                                #print "iiiiiiinnnnn", loopEntrie, loopEntrie[4].upper()
                                #print sortedSingleLoopDict[index]
                                if not sortedSingleLoopDict[index][8]:
                                    # print "1", sortedSingleLoopDict[index][8], sortedSingleLoopDict[index][2]
                                    sortedSingleLoopDict[index][2] *= self.GPCRscore
                                    # print "2", sortedSingleLoopDict[index][8], sortedSingleLoopDict[index][2]
                                sortedSingleLoopDict[index].append(True)
                                    # print sortedSingleLoopDict[index][2]
                            else :
                                sortedSingleLoopDict[index].append(False)

                        #print sortedSingleLoopDict

                        sortedSingleLoopDict = filter(lambda loop: loop[5] < self.numclashes, sortedSingleLoopDict)
                        sortedSingleLoopDict.sort(key=lambda loop: (loop[2], loop[5]), cmp =self.compareLoops )
                        sortedSingleLoopDict.reverse()

                        #print sortedSingleLoopDict


                        for resultIndex in range(min(20, len(sortedSingleLoopDict))) :
                            result = sortedSingleLoopDict[resultIndex]
                            #print result
                            # [0] Erweiterung, Loop, Loopauswahl, Datenbank [1] Score, [2] Clashes, [3] Sequenz, [4] Template, [5] Position(Datenbank), [6] GPCR True/False, [7] SequenzIdentitaet, [8] Spezies True/False, [9] memDB True/False
                            loop.append(["%i,%i;%i;%i;%i" % (i,i,j,result[0],memDB), result[2], result[5], result[3], result[4], result[6], result[9], round(result[7],2), result[8], memDB])
                            # print loop
                        temp_loopDict[j] = loop
                        #temp_loopDict[j] = filter(lambda loop: not (loop[0].split(',')[0] == 0 and loop[2] != 0), temp_loopDict[j])

            return temp_loopDict, singleLoopDict

        # print loopDict
        outputDir1 = self.output_dir + "membranDB"
        loopDict1, singleLoopDict1 = firstsort( outputDir1, memDB=False )
        outputDir2 = self.output_dir + "GPCRDB"
        loopDict2, singleLoopDict2 = firstsort( outputDir2, memDB=True)

        if not singleLoopDict1 == {} or not singleLoopDict2 == {} :


            # print singleLoopDict1
            # print '===================================================================='
            # print singleLoopDict2
            # print singleLoopDict
            loopDict = {}
            if singleLoopDict1 == {}:
                loopDict = loopDict2
            elif singleLoopDict2 == {}:
                loopDict = loopDict1
            else:
                for modelnum in [0,1,2,3,4,5,6]:#loopDict2:
                    # ,4,5
                    # print '================================================================'
                    # print modelnum
                    # falls einer der Loops (ICL1) in einem der dicts nicht vorhanden ist, abfrage mit get und []
                    #loopDict[modelnum] = loopDict1 + loopDict2
                    loopDict[modelnum] = loopDict1.get(modelnum, []) + loopDict2.get(modelnum, [])
            #print loopDict
            #Ergebnis wird in BestResultsDict.json ausgegeben
            with open(os.path.join(self.output_dir, "BestResultsDict.json"), 'w') as fp :
                json.dump(loopDict, fp)

            # wird spaeter sortiert fuer einzelne pdb-files in pdbLoop()
            #print loopDict
            #print '======================='
            self.sortedPdbLoopDict = loopDict.copy()
            # print loopDict

            #topTenResult entaelt die Loopinformationen der 10 besten gesamt Ergebnisse
            # --> benutzt fuer zusammenbauen der 10 besten pdb-files
            self.topTenResult = {}

            #sortiert nach hoestem Score die Liste der 3 besten berechneten Loops pro Loop (12 Ergebnisse pro Loop) und dreht sie um, d.h bester Score zu erst

            for key in loopDict :
                loopDict[key] = filter(lambda loop: loop[2] < self.numclashes, loopDict[key])
                loopDict[key].sort(key=lambda loop: (loop[1], loop[2]), cmp =self.compareLoops )
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

            self.avgScoreList = []
            scoreList = []

            for key in self.topTenResult :
                for loop in self.topTenResult[key] :
                    scoreList.append(loop[1])
                avgScore = round(sum(scoreList)/len(scoreList),3)
                self.avgScoreList.append((key, avgScore))

            with open(os.path.join(self.output_dir, "AvgScore.json"), 'w') as fp :
                json.dump(self.avgScoreList, fp)

            # self.insertLoop()

            self.pdbLoop()

        else :
            with open ('FailedLoop.txt' , 'w') as fp :
                fp.write( 'Es konneten keine Loops berechnet werden')

    def pdbLoop( self ):

        # print '---------------'
        # print self.sortedPdbLoopDict
        #shutil.copyfile('../../../index.html', os.path.join(self.output_dir, "index.html" ))

        relResultLoopsList = []
        allResultLoopList = ['ICL1', 'ICL2', 'ICL3', 'ECL1', 'ECL2', 'ECL3' , '8']

        templateHtmlString = ''

        if self.sortedPdbLoopDict == {} :
            return

        # print self.sortedPdbLoopDict
        # print '==========================='



        n=0
        k=0

        #loecht eintraege je looptyp, wenn tamplate doppelt vorkommt (z.b.4z36)
        #print self.sortedPdbLoopDict
        for key in self.sortedPdbLoopDict :
            noDoubelList = []
            templateList = []
            # print loopElements
            for loopElement in self.sortedPdbLoopDict[key] :
                extension, loopName, numLoop, numDB = loopElement[0].split(';')
                extensionSingle1, extensionSingle2 = extension.split(',')
                extensionSingle1 = int(extensionSingle1)
                posNum, posLet = loopElement[5].split(':')
                posNum1, posNum2 = filter(None, posNum.split('-'))
                posNum1 = int(posNum1)
                # print extensionSingle1
                # print posNum1
                seqPos = posNum1+extensionSingle1
                # print seqPos
                # print "--------"
                if [loopElement[4],seqPos] in templateList :
                    n = n+1
                    #print (loopElement[4],seqPos)
                else :
                    templateList.append([loopElement[4], seqPos])
                    noDoubelList.append(loopElement)
                    n = n+1
            self.sortedPdbLoopDict[key] = noDoubelList





        # print self.sortedPdbLoopDict
        # print n
        # print k


        for key in self.sortedPdbLoopDict:
            self.sortedPdbLoopDict[key] = filter(lambda loop: loop[2] < self.numclashes, self.sortedPdbLoopDict[key])
            self.sortedPdbLoopDict[key].sort(key=lambda loop: (loop[1], loop[2]), cmp =self.compareLoops )
            self.sortedPdbLoopDict[key].reverse()
            self.sortedPdbLoopDict[key] = self.sortedPdbLoopDict[key][0:5] # er werden die ersten 5 genommen

        with open(os.path.join(self.output_dir, "SortedPdBLoopDict.json"), 'w') as fp :
            json.dump(self.sortedPdbLoopDict, fp)

        resultLoopNames = []
        oriSequenceList = []
        makeLoopList = []
        k = 0
        for key in self.sortedPdbLoopDict :
            loopLineList = []
            loopStartEndList = []
            for loop in self.sortedPdbLoopDict[key] :
                extension, loopName, numLoop, numDB = loop[0].split(';')

                loopName = int(loopName)
                numLoop = int(numLoop)
                numDB = int(numDB)
                #da extension ein Doupel ist (z.B. 0,0), wird es aufgeteilt
                extensionSingle1, extensionSingle2 = extension.split(',')
                extensionSingle1 = int(extensionSingle1)
                extensionSingle2 = int(extensionSingle2)

                loopFile = []
                if numDB == 1:
                    databaseOri = "GPCRDB"
                else:
                    databaseOri = "membranDB"
                with open(self.subdir(databaseOri+"/link_it_%i,%i/link_it_%i" % (extensionSingle1,extensionSingle2,loopName)) + "/" + os.path.splitext(os.path.basename(self.pdb_file))[0] + "_input_out_linker2.pdb", 'r') as fp :
                    loopFile = fp.readlines()

                modelBegin = -1
                modelEnd = -1
                for i, line in enumerate(loopFile) :
                    if line[0:5] == 'MODEL' and line.split()[1] == str(numLoop) :
                        modelBegin = i
                    if line[0:6] == 'ENDMDL' and modelBegin != -1 :
                        modelEnd = i
                        break

                loopLineList.append(loopFile[modelBegin + 5 : modelEnd - 4])
                loopStart = int(loopFile[modelBegin + 5][22:26])
                loopEnd = int(loopFile[modelEnd - 4][22:26])
                loopStartEndList.append((loopStart, loopEnd))
                # print modelBegin
                # print modelEnd
                # print loopLineList
                oriModel = loopStart + extensionSingle1



                # fuer Name vom  Loop
                with open( self.loop_file, 'r' ) as oriLoop:
                    for line in oriLoop:
                        if line[0] == '#':
                            continue
                        match = re.search("(\w*):\s*start=(\d*):\s*end=\d*:\s*(\w*)",line)
                        if match:
                            oriloopName, start, oriSequence= match.groups()
                            #print oriloopName
                            #print start
                            if int(start) == oriModel :
                                resultLoopNames.append(oriloopName)
                                oriSequenceList.append(oriSequence)



            #print resultLoopNames
            #print oriSequenceList



            # for n, i in enumerate(resultLoopNames) :
            #     if i == '8' :
            #         resultLoopNames[n] == 'TM7_H8'

            resultLoopNames = [ 'TM7_H8' if x == '8' else x for x in resultLoopNames ]


            #print relResultLoopsList



            with open ( self.ori_file, 'r' ) as oriFile:
                oriFileContent = oriFile.readlines()

            for j, result in enumerate(loopLineList) :
                new_pdb_file =  (resultLoopNames[k] + '_%i.pdb' %j)
                makeLoopList.append(resultLoopNames[k] + '_%i.pdb' %j)


                #new_pdb_file = (os.path.splitext(self.ori_file)[0] + '_out_%i.pdb' %j)
                with open (new_pdb_file, 'w') as fp1 :
                    loopStart = loopStartEndList[j][0]
                    loopEnd = loopStartEndList[j][1]
                    for line in [l for l in oriFileContent if l[0:4] == 'ATOM' and loopStart - 3 <= int(l[22:26])  < loopStart] :
                        fp1.write(line)
                    for line in result :
                        fp1.write(line)
                    for line in [l for l in oriFileContent if l[0:4] == 'ATOM' and loopEnd  <= int(l[22:26])  < loopEnd + 3] :
                        fp1.write(line)

                k += 1
                #print k

        # print resultLoopNames
        # print allResultLoopList

        with open( self.loop_file, 'r' ) as oriLoop:
                for line in oriLoop:
                    if line[0] == '>' :
                        relResultLoopsList.append(line[1:-1])

        for name in allResultLoopList :
                if name in resultLoopNames :
                    relResultLoopsList.append((name,True))
                else :
                    relResultLoopsList.append((name,False))

        # print relResultLoopsList

        relResultLoopsList = [ ('TM7_H8',x[1]) if x[0] == '8' else x for x in relResultLoopsList ]

        # print '================='
        # print relResultLoopsList

        with open(os.path.join(self.output_dir, "RelResultLoopsList.json"), 'w') as fp :
                json.dump(relResultLoopsList, fp)

        #ErgebnissTabelle erstellen

        i = 0
        tableList = self.sortedPdbLoopDict.copy()
        # print self.sortedPdbLoopDict
        nameList = []

        lastName = ''
        for name in resultLoopNames :
            if name == lastName :
                continue
            nameList.append(name)
            lastName = name

        # for loops in tableList.values() :
        #     for loop in loops :
        #         loop.insert(4, oriSequenceList[i])
        #         i+=1

        # print tableList
        # print '==================================================='
        # print self.oriSeqDict

        for key, loops in tableList.iteritems() :
            for loop in loops :
                loop.insert(4,self.oriSeqDict[int(loop[0][0][0])][key])

        # print '===================================================='
        # print tableList
   #     for i in range(7):
  #          if i in tableList:
 #               fixedTableList[]


        with open ( 'Result_Table.csv', 'w' ) as fp1 :
            fp1.write('Loop ,# , Ext, GPCR, Species, Score, Sequence, Templ seq, Seq ident, Clashes, PDB-ID, Templ pos\n' )
            j = 0
            for name in (nameList) :
                fp1.write(name)
                #print tableList
                while j not in tableList or tableList[j] == []:
                    j += 1
                for i,  loop in enumerate(tableList[j]) :
                    # print loop
                    # score erhoeungen zurueck aendern
                    newscore = loop[1]
                    if loop[9]:#Species
                        newscore = newscore/self.Speciesscore
                    elif loop[7]:#GPCR
                        newscore = newscore/self.GPCRscore

                    fp1.write(','+ ('%i' %i) + ',' + loop[0][0] + ',' + str(loop[7]) +  ',' + str(loop[9]) + ',' + str(newscore) + ',' + loop[4] + ',' + loop[3] + ',' + str(loop[8]) + ',' + str(loop[2]) + ',' + loop[5] + ',' + loop[6] + '\n')
                j += 1

        with open ( 'Result.txt', 'w') as fp2 :
            j = 0
            for name in (nameList) :
                while j not in tableList or tableList[j] == [] :
                    j += 1
                for i,  loop in enumerate(tableList[j]) :
                    # print loop
                    newscore = loop[1]
                    if loop[9]:#Species
                        newscore = newscore/self.Speciesscore
                    elif loop[7]:#GPCR
                        newscore = newscore/self.GPCRscore
                    fp2.write( name + ',' + name + ('_%i' %i) + ',' + loop[0][0] + ',' + str(loop[7]) + ',' + str(loop[9]) + ',' + str(newscore) + ',' + loop[4] + ',' + loop[3] + ',' + str(loop[8]) + ',' + str(loop[2]) + ',' + loop[5] + ',' + loop[6] + ';')
                j += 1

        result = ""
        j = 0
        for name in (nameList) :
            while j not in tableList or tableList[j] == [] :
                j += 1
            for i,  loop in enumerate(tableList[j]) :
                # print loop
                newscore = loop[1]
                if loop[9]:#Species
                    newscore = newscore/self.Speciesscore
                elif loop[7]:#GPCR
                    newscore = newscore/self.GPCRscore
                result +=  name + ',' + name + ('_%i' %i) + ',' + loop[0][0] + ',' + str(loop[7]) + ',' + str(loop[9]) + ',' + str(newscore) + ',' + loop[4] + ',' + loop[3] + ',' + str(loop[8]) + ',' + str(loop[2]) + ',' + loop[5] + ',' + loop[6] + ';'
            j += 1

        with open ( 'Result_Tabel.txt', 'w') as fp2 :
            j = 0
            for name in (nameList) :
                while j not in tableList or tableList[j] == [] :
                    j += 1
                for i,  loop in enumerate(tableList[j]) :
                    # print loop
                    newscore = loop[1]
                    if loop[9]:#Species
                        newscore = newscore/self.Speciesscore
                    elif loop[7]:#GPCR
                        newscore = newscore/self.GPCRscore
                    fp2.write( name + ',' + name + ('_%i' %i) + ',' + loop[0][0] + ',' + str(loop[7]) + ',' + str(loop[9]) + ',' + str(newscore) + ',' + loop[4] + ',' + loop[3] + ',' + str(loop[8]) + ',' + str(loop[2]) + ',' + loop[5] + ',' + loop[6] + ';')
                j += 1

        resultTabel = "Loop ,# , Ext, GPCR, Species, Score, Sequence, Templ seq, Seq ident, Clashes, PDB-ID, Templ pos,;"
        j = 0
        PDBaddress = 'http://www.rcsb.org/pdb/explore/explore.do?structureId='
        for name in (nameList) :
            while j not in tableList or tableList[j] == [] :
                j += 1
            for i,  loop in enumerate(tableList[j]) :
                # print loop
                newscore = loop[1]
                if loop[9]:#Species
                    newscore = newscore/self.Speciesscore
                elif loop[7]:#GPCR
                    newscore = newscore/self.GPCRscore
                resultTabel += name + ',' + ('%i' %i) + ',' + loop[0][0] + ',' + str(loop[7]) + ',' + str(loop[9]) + ',' + str(newscore) + ',' + loop[4] + ',' + loop[3] + ',' + str(loop[8]) + ',' + str(loop[2]) + ',' + ('<a target="_blank" href="' + PDBaddress + loop[5] + '">' + loop[5] + '</a>') + ',' + loop[6] + ';'
            j += 1



        # print makeLoopList

        shutil.copyfile(downloader_js_file, os.path.join(self.output_dir, "downloader.js" ))
        #html.index erstellen
        shutil.copyfile(ngl_js_file, os.path.join(self.output_dir, "ngl.js" ))
        # tmpl_name="ngl.embedded.min.js"
        # tmpl_file = os.path.join( tmpl_dir, tmpl_name )
        new_pdb_file = (self.ori_file)
        templateHtmlString += '"' + os.path.basename(new_pdb_file) + '"'
        #print templateHtmlString

        templateHtml = ''
        with open (index_html_file, 'r') as fp :
            templateHtml = fp.read()

        loopCandidate = ["ECL1_0","ECL1_1","ECL1_2","ECL1_3","ECL1_4","ECL2_0","ECL2_1","ECL2_2","ECL2_3","ECL2_4",
                         "ECL3_0","ECL3_1","ECL3_2","ECL3_3","ECL3_4",
                         "ICL1_0","ICL1_1","ICL1_2","ICL1_3","ICL1_4","ICL2_0","ICL2_1","ICL2_2","ICL2_3","ICL2_4",
                         "ICL3_0","ICL3_1","ICL3_2","ICL3_3","ICL3_4",
                         "TM7_H8_0","TM7_H8_1","TM7_H8_2","TM7_H8_3","TM7_H8_4"]

        posHelix = []
        fileContent = []
        with open( self.ori_file, 'r') as fp :
            fileContent = fp.readlines()

        for u, line in enumerate(fileContent) :
            if line[0:4] == 'ATOM' and fileContent[u+1][0:4] == 'ATOM':
                startHelix = int(line[22:26])
                posHelix.append(startHelix)
                break
            else :
                continue

        for w, line in enumerate(fileContent) :
            if line[0:4] == 'ATOM' :
                startHelix = int(line[22:26])
                if fileContent[w+1][0:4] == 'ATOM' :
                    nextHelix = int(fileContent[w + 1][22:26])
                    if nextHelix == startHelix or nextHelix == startHelix + 1:
                        continue
                    else :
                        posHelix.append(startHelix)
                        posHelix.append(nextHelix)

        for t, line in enumerate(fileContent)  :
            if line[0:4] == 'ATOM' and fileContent[t+1][0:4] != 'ATOM':
                startHelix = int(line[22:26])
                posHelix.append(startHelix)
            else :
                continue


        # print posHelix

        loopPosList = ["loopPos1","loopPos2","loopPos3","loopPos4","loopPos5","loopPos6","loopPos7"]
        templateDict = {}

        for v in loopCandidate :
            if (v + '.pdb') in makeLoopList :
                templateDict[v] = '"'+ str(v + '.pdb')+'"'
            else :
                templateDict[v] = ""

        # print "================="
        # print loopPosList
        # speichert positionen der helix fuer die html, nimmt dafuer die luecken im pdb.file
        q = 0
        for p in loopPosList :
            templateDict[p] = '"' + str(posHelix[q]) + '-' + str(posHelix[q+1]) + '"'
            q += 2
        # print "====================="
        # print templateDict

        templateDict["namepdbfiles"] = templateHtmlString
        templateDict["namepdbfiles2"] = templateHtmlString
        templateDict["result"] = '"' + result + '"'
        templateDict["resultTabel"] = "'" + resultTabel + "'"
        #print templateDict

        # print self.loopPosList
        templateHtml = templateHtml.format(**templateDict)

        #print templateHtml


        with open (os.path.join(self.output_dir, "index.html" ), 'w') as fp :
            fp.write(templateHtml)


    #fuegt Loops in pdb-file ein
    #--> Ergebnis sind die 10 besten pdb-files
    def insertLoop( self ):

        #shutil.copyfile('../../../index.html', os.path.join(self.output_dir, "index.html" ))

        templateHtmlString = ''

        if self.topTenResult == {} :
            return

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

                loopLineList.append(loopFile[modelBegin + 5 : modelEnd - 4])

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
            shutil.copyfile('../../../ngl.js', os.path.join(self.output_dir, "ngl.js" ))
            # tmpl_name="ngl.embedded.min.js"
            # tmpl_file = os.path.join( tmpl_dir, tmpl_name )
            templateHtmlString += '"' + os.path.basename(new_pdb_file) + '",\n'
            # print templateHtmlString

            templateHtml = ''
            with open ('../../../index.html', 'r') as fp :
                templateHtml = fp.read()

            #print self.loopPosList
            templateHtml = templateHtml.format(namepdbfiles = templateHtmlString,
                                               loopPos1 = '"'+ str(self.loopPosList[0][0]) + '-' + str(self.loopPosList[0][1])+'"',
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
        _( "names", type="list", nargs="*", default=None ),
        _( "gpcrDB", type="boolean", default=False ),
        _( "last_res", type="int", default=9999 ),
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
            # print LINKIT_CMD
            # print '============'
            # print LINKIT_DIR2
            # print 'res1:  ' , res1
            # print 'res2:  ' , res2
            # print 'seq:  ' , seq
            # print 'last:  ' , self.last_res
            if int(res1.split(":")[0]) <= self.last_res and int(res2.split(":")[0]) <= self.last_res :
                link_it = LinkIt(
                    self.pdb_file, res1, res2, seq, memdb=not self.gpcrDB, GPCRdb=self.gpcrDB,
                    **copy_dict( kwargs, run=False, debug=True,
                                 output_dir=self.subdir("link_it_%i" % i) )
                )
                self.output_files += link_it.output_files
                self.sub_tool_list.append( link_it )
                self.link_it_list.append( link_it )
                #print self.link_it_list


    def func( self ):
        self.log( "%i linkit runs" % len( self.input ) )
        for i, link_it in enumerate ( self.link_it_list ):
            # print 'res1:  ' , link_it.res1
            # print 'res2:  ' , link_it.res2
            # print 'seq:  ' , link_it.seq
            try:
                link_it()
            except:
                pass
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
        _( "pdb_file", type="file", ext="pdb", label="PDB file",
            help="The input structure." ),
        _( "mrc_file", type="file", ext="mrc", label="Density file",
            help="The input density." ),
        _( "res1", type="sele", label="Stem residue 1",
            help="N-terminal stem residue and chain, '123:A'." ),
        _( "res2", type="sele", label="Stem residue 2",
            help="C-terminal stem residue." ),
        _( "seq", type="str", label="Fragment sequence",
            help="One-letter code of the linker amino acids." ),
        _( "resolution", type="float", fixed=True , range=[0.5,20], default=5,
            precision=1, label="Map resolution",
            help="Used for filtering the map." ),
        _( "memdb", type="bool", label="MembraneDB",
            help="Show only results from membrane proteins.", default=False ),
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
        self.res1['chain'] = self.res1['chain'].upper()
        self.res2['chain'] = self.res2['chain'].upper()
        if self.res1['resno'] > self.res2['resno']:
            self.res1, self.res2 = self.res2, self.res1
        self.link_it = LinkIt(
            self.edited_pdb_file, self.res1, self.res2, self.seq, self.memdb,
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
            #print fn
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
