from __future__ import with_statement
from __future__ import division

import os
import itertools
import json
import collections
import xml.etree.ElementTree as ET

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
        _( "seq", type="str", label="Sequence",
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

        oripdb=npdb.sele(chain=self.res1['chain'],resno=[self.res1['resno']-2,self.res2['resno']+2] ,invert=True)
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


class MultiLinkIt( PyTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "input", type="list", nargs=3, action="append",
           help="sele,sele,str" ),
        _( "names", type="list", nargs="*", default=None )
    ]
    out = []
    tmpl_dir = TMPL_DIR
    provi_tmpl = "multi_link_it.provi"

    def _init( self, *args, **kwargs ):
        self.link_it_list = []
        for i, linker_args in enumerate( self.input ):
            res1, res2, seq = linker_args
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
        for link_it in self.link_it_list:
            link_it()

    def _post_exec( self ):
        linker_list = []
        for i, link_it in enumerate( self.link_it_list ):
            linker_list.append({
                "i": i,
                "name": self.names[i] if self.names else "Linker",
                "json_file": self.relpath( link_it.json_file ),
                "pdb_linker_file3": self.relpath( link_it.pdb_linker_file3 )
            })
        self._make_provi_file(
            use_jinja2=True,
            pdb_file=self.relpath( self.pdb_file ),
            linker_list=linker_list
        )


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
        print len(self.seq)
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
