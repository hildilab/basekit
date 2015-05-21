from __future__ import with_statement
from __future__ import division
                                                                                                                                                                                                                                                                                                                                                                                                                        

import os
import itertools
import operator
import json
import collections
import xml.etree.ElementTree as ET
import timeit
import utils.path
from utils.timer import Timer
from utils import copy_dict, iter_stride, dir_walker
from utils.tool import _, _dir_init, CmdTool, PyTool, ProviMixin
from utils.numpdb import NumPdb, numsele
from utils.mrc import get_mrc_header, getMrc
import numpy as np
import time
import provi_prep as provi
from spider import LoopCrosscorrel
from pdb import PdbEdit, SplitPdbSSE, LoopDelete ,PdbSplit,get_tree
import utils.numpdb as numpdb

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "linker" ) 

def LINKIT_DIR():
    return os.environ.get("LINKIT_DIR", "")

def LINKIT_CMD():
    return os.path.join( LINKIT_DIR(), "Link_It_dos2n.exe" )



class LinkIt( CmdTool, ProviMixin ):
    """Please upload the PDB file, define the N- and C-terminal stem-residues (e.g.: 10:A, 16:A) and provide the missing sequence in 1-letter code (ACDEF)."""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "res1", type="sele",help="resno:chain, i.e. 10:A" ),
        _( "res2", type="sele" ),
        _( "seq", type="str" ),
        _( "max_loops", type="int", range=[0, 500], default=500 )
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
        if self.res1['resno'] > self.res2['resno']:
            self.res1, self.res2 = self.res2, self.res1
        self.cmd = [ "wine", LINKIT_CMD(), self.kos_file, self.bin_file, "t" ]

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
                    fp_out.write( line )
                if line.startswith("ATOM"):
                    line = line[0:6] + ( "% 5i" % atom_i ) + line[11:]
                    if line[22] == "X":
                        if not stems:
                            continue
                        if line[12:16] not in backbone:
                            continue
                        if line[24] == " ":
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
                    fp_out.write( line )

    def _split_loop_file( self ):
        PdbSplit(
            self.pdb_linker_file2, output_dir=self.loop_dir,
            backbone_only=True, resno_ignore=[ 1000, 2000 ], zfill=3,max_models=self.max_loops
        )

    def _find_clashes( self ):
        npdb = NumPdb( self.pdb_file, { "backbone_only": True } )
        protein_tree = get_tree(npdb['xyz'])
        atom_resno_list = npdb.get('resno')

        model_clash_count = {}

        for m, file_path in dir_walker( self.loop_dir, ".*\.pdb" ):
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
                if e not in (self.res1['resno'], self.res2['resno']):
                    clashes.append(e)

            model_no = int(
                utils.path.basename( file_path ).split('_')[0]
            )
            model_clash_count[ model_no ] = len( clashes )

        return model_clash_count

    def _make_linker_json( self, compact=False ):
        linker_dict = {}
        model_clash_count = self._find_clashes()

        with open( self.txt_file, "r" ) as fp:
            fp.next()
            fp.next()
            for i, d in enumerate( iter_stride( fp, 4 ), start=1 ):
                if i<=self.max_loops:
                    linker_dict[ i ] = [
                        float(d[0]), float(d[1]),
                        str(d[2].strip()),
                        str(d[3].strip()),
                        model_clash_count[ i ]
                    ]

        with open( self.json_file, "w" ) as fp:
            if compact:
                json.dump( linker_dict, fp, separators=(',', ':') )
            else:
                json.dump( linker_dict, fp, indent=4 )

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
           self.link_it_list.append( link_it )
   def func( self ):
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


class LinkItDensity( PyTool, ProviMixin ):
    """Please choose a cutoff  and denote the resolution of your map. Then define the stem-residues (e.g.: 10:A, 16:A) and provide the missing sequence in 1-letter code (ACDEF). """
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "mrc_file", type="file", ext="mrc" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "seq", type="str" ),
        _( "resolution", type="float", fixed=True ),
        _( "cutoff", type="float" ),
        _( "max_loops", type="int", range=[0, 500], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it_density.provi"
    #start=timeit.timeit()
    def _init( self, *args, **kwargs ):
        if self.res1['resno'] > self.res2['resno']:
            self.res1, self.res2 = self.res2, self.res1
        self.link_it = LinkIt(
            self.edited_pdb_file, self.res1, self.res2, self.seq,
            **copy_dict( kwargs, run=True, output_dir=self.subdir("link_it") )
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
        print shx
        PdbEdit(
            self.pdb_file, shift=[ shx, shy, shz ]
        )

        self.link_it()
        self.loop_correl()
        print len(self.seq)
    def _post_exec( self ):
        self._make_correl_json()
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            # pdb_file=self.relpath( self.loop_correl.spider_shift.edited_pdb_file ),
            mrc_file=self.relpath( self.mrc_file ),
            # mrc_file=self.relpath( self.loop_correl.spider_shift.map_shift ),
            cutoff=self.cutoff,
            # box_mrc_file=self.relpath(
            #     self.loop_correl.spider_reconvert.mrc_file
            # ),
            # box_ori_mrc_file=self.relpath(
            #     self.loop_correl.spider_reconvert.mrc_ori_file
            # ),
            pdb_linker_file3=self.relpath(
                self.loop_correl.ori_pdb_linker_file3 ),
            linker_correl_file=self.relpath( self.linker_correl_file )
        )
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
        linker_correl_dict = {}
        for k, v in cc_dict.iteritems():
            linker_correl_dict[k] = [ v ] + li_dict[k]
        with open( self.linker_correl_file, "w" ) as fp:
            if compact:
                json.dump( linker_correl_dict, fp, separators=(',', ':') )
            else:
                json.dump( linker_correl_dict, fp, indent=4 )
    
    #end = timeit.timeit()
    #print "zEdIT", end - start
class LnkItVali(PyTool, ProviMixin):
    args = [
        _( "dataset_dir", type="dir" ),
        _( "cutoff", type="float", default=5 ),
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
    _( "cutoff", type="float" )
        ]
    def func( self, *args, **kwargs ):
        
        npdb=numpdb.NumPdb( self.pdb_file)

        for z, numa in enumerate (npdb.iter_chain()):
            cha=numa.get('chain')[0]
            rel=numa.get('resno')[-1]
            print rel
            print cha
            os.mkdir(cha)
            for x in range (28,36,1):
                ch=numa.get('chain')[0]
                gum= "%s/%s" % (cha,x)
                os.mkdir(gum)
                for i,numa in enumerate (npdb.iter_resno(chain=cha)):
                   
                    
                   
                    res1=numa.get('resno')[0]
                    if res1 % 5 == 0 :#-->jeder 7. loop
                        res2=res1+x-1
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
                           # start = timeit.timeit()
                            #print 'richtige?',numa.sequence()
                            
                            LinkItDensity (name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=di)
                            #except:
                             #   print 'linikt error'
                            #    neuen loop suchen
            #break
                            
class CutPDB2 (PyTool):
    args = [
    _( "pdb_file", type="file"),
    _( "mrc_file", type="file"),
    _( "resolution", type="float", range=[1, 10], step=0.1 ),
    _( "cutoff", type="float" )
        ]
    def func( self, *args, **kwargs ):
        
        npdb=numpdb.NumPdb( self.pdb_file)

        for z, numa in enumerate (npdb.iter_chain()):
            #chain=numa.get('chain')[0]
            #print chain
            #if chain orinpdb=numpdb.NumPdb(ori, {
                print 'numa', numa
                cha=numa.get('chain')[0]
                rel=numa.get('resno')[-1]
                print rel
                print cha
                os.mkdir(cha)
                for x in range (5,35,1):
                    ch=numa.get('chain')[0]
                    gum= "%s/%s" % (cha,x)
                    os.mkdir(gum)
                    orinpdb=numpdb.NumPdb(self.pdb_file, {'chain':cha})
                    for i,numb in enumerate (orinpdb.iter_resno(chain=cha)):
                       
                        
                       
                        res1=numb.get('resno')[0]
                        if res1 % 7 == 0 :#-->jeder 7. loop
                            res2=res1+x-1
                            print res1
                            test=npdb.copy(resno=[res1,res2],chain=cha)
                            oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=True)
                            if rel-x+1>=res1:
                                di="%s/%s_%s" % (gum,res1,res2)
                                os.mkdir(di)
                                name= "%s/%s_%s.%s" % (di,res1,res2,'pdb')
                                name2= "%s/%s_%s_%s.%s" % (di,'ganz',res1,res2,'pdb')
                                print name2
                                #test.write(name)
                                #npdb.write(name2,sele=oripdb)
                                seq=test.sequence()
                                ires1="%s:%s" % (res1-1,cha)
                                ires2="%s:%s" % (res2+1,cha)
                                print ires1, ires2
                                print 'sequence',seq
                               # start = timeit.timeit()
                                #print 'richtige?',numa.sequence()
                                try:
                                    LinkItDensity (name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=di)
                                except:
                                    print 'linikt error'
                                #    neuen loop suchen
                #break
                                

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
                    
#class Timer:    
#    def __enter__(self):
#        self.start = time.clock()
#        return self
#
#    def __exit__(self, *args):
#        self.end = time.clock()
#        self.interval = self.end - self.start
#/home/jochen/work/fragfit/validation/dataset/EMD-5249/pieces/3IZM_B_705-708/loop_correl/crosscorrelation
class Cut_analyse_PDB2 (PyTool):
    args = [
    _( "pdb_file", type="file"),
    _( "mrc_file", type="file"),
    _( "resolution", type="float", range=[1, 10], step=0.1 ),
    _( "cutoff", type="float" )
        ]
    def func( self, *args, **kwargs ):
        npdb=numpdb.NumPdb( self.pdb_file)
        cha='R' 
        rel=198
        os.mkdir(cha)
        for x in range (5,36,1):
            ch='R'
            gum= "%s/%s" % (cha,x)
            os.mkdir(gum)
            res1=198
            res2=res1+x-1
            test=npdb.copy(resno=[res1,res2],chain=cha)
            oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=True)
            di="%s/%s_%s" % (gum,res1,res2)
            os.mkdir(di)
            name= "%s/%s_%s.%s" % (di,res1,res2,'pdb')
            name2= "%s/%s_%s_%s.%s" % (di,'ganz',res1,res2,'pdb')
            test.write(name)
            npdb.write(name2,sele=oripdb)
            seq=test.sequence()
            ires1="%s:%s" % (res1-1,cha)
            ires2="%s:%s" % (res2+1,cha)
            with Timer() as t:
                print len(seq)
                try:
                   # LinkIt_analyse (name2,ires1,ires2,seq)
                #print name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,di
                    LinkItDensity2_analyse (name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=di)
            #print t
                except:
                    print 'linikt error'
            #    neuen loop suchen
            #
                
            #except:
            #    print 'aha'
            #    continue
            #end = timeit.timeit()
            #print end - start
            #print t.interval
        #else:
            #    continue
                #        
class LinkItDensity2_analyse( PyTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "mrc_file", type="file", ext="mrc" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "seq", type="str" ),
        _( "resolution", type="float", fixed=True ),
        _( "cutoff", type="float" ),
        _( "max_loops", type="int", range=[0, 500], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it_density.provi"
    #start=timeit.timeit()
    def _init( self, *args, **kwargs ):
        if self.res1['resno'] > self.res2['resno']:
            self.res1, self.res2 = self.res2, self.res1
        self.link_it = LinkIt_analyse(
            self.edited_pdb_file, self.res1, self.res2, self.seq,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("link_it") )
        )
        self.loop_correl = LoopCrosscorrel2(
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

    def func( self ):
        boxsize = getMrc(self.mrc_file, 'nx' )
        originx = abs(getMrc(self.mrc_file, 'nxstart' ))
        originy = abs(getMrc(self.mrc_file, 'nystart' ))
        originz = abs(getMrc(self.mrc_file, 'nzstart' ))
        size = getMrc(self.mrc_file, 'xlen' )
        pixelsize = size / boxsize
        shx = (originx - (boxsize / 2)) * pixelsize
        shy = (originy - (boxsize / 2)) * pixelsize
        shz = (originz - (boxsize / 2)) * pixelsize
        PdbEdit(
            self.pdb_file, shift=[ shx, shy, shz ]
        )

        self.link_it()
        self.loop_correl()
        print len(self.seq)
    def _post_exec( self ):
        self._make_correl_json()
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            # pdb_file=self.relpath( self.loop_correl.spider_shift.edited_pdb_file ),
            mrc_file=self.relpath( self.mrc_file ),
            # mrc_file=self.relpath( self.loop_correl.spider_shift.map_shift ),
            cutoff=self.cutoff,
            # box_mrc_file=self.relpath(
            #     self.loop_correl.spider_reconvert.mrc_file
            # ),
            # box_ori_mrc_file=self.relpath(
            #     self.loop_correl.spider_reconvert.mrc_ori_file
            # ),
            pdb_linker_file3=self.relpath(
                self.loop_correl.ori_pdb_linker_file3 ),
            linker_correl_file=self.relpath( self.linker_correl_file )
        )
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
        linker_correl_dict = {}
        for k, v in cc_dict.iteritems():
            linker_correl_dict[k] = [ v ] + li_dict[k]
        with open( self.linker_correl_file, "w" ) as fp:
            if compact:
                json.dump( linker_correl_dict, fp, separators=(',', ':') )
            else:
                json.dump( linker_correl_dict, fp, indent=4 )
    
    #end = timeit.timeit()
    #print "zEdIT", end - start

class CutPDB3 (PyTool):
    args = [
    _( "pdb_file", type="file"),
    _( "mrc_file", type="file"),
    _( "resolution", type="float", range=[1, 10], step=0.1 ),
    _( "cutoff", type="float" )
        ]
    def func( self, *args, **kwargs ):
        
        npdb=numpdb.NumPdb( self.pdb_file)

        for z, numa in enumerate (npdb.iter_chain()):
            cha=numa.get('chain')[0]
            rel=numa.get('resno')[-1]
            print rel
            print cha
            #if cha not in ('1','0'):
            try:
                os.mkdir(cha)
            except:
                continue
            for x in range (7,36,1):
                ch=numa.get('chain')[0]
                gum= "%s/%s" % (cha,x)
                
                os.mkdir(gum)
                if ch not in ('A','0','1','2','3','4','5','6','7','8','9'):
                    
                    for i,numa in enumerate (npdb.iter_resno(chain=cha)):
                        #numa.next()
                        #numa.next()
                        #numa.next()
                        res1=numa.get('resno')[0]
                        if res1 % 7 == 0 :#-->jeder 7. loop
                            res2=res1+x-1
                            print res1
                            test=npdb.copy(resno=[res1,res2],chain=cha)
                            oripdb=npdb.sele(resno=[res1,res2],chain=cha, invert=True)
                            if rel-x+1>=res1:
                                di="%s/%s_%s" % (gum,res1,res2)
                                try:
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
                                        LinkItDensity (name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=di)
                                    except:
                                        print 'linikt error'
                                except:
                                    continue
                    else:
                        continue
                                #    neuen loop suchen
                        
class CalcSheets (PyTool):
        args = [
    _( "pdb_file", type="file"),
    _( "mrc_file", type="file"),
    _( "resolution", type="float", range=[1, 10], step=0.1 ),
    _( "cutoff", type="float" )
        ]
        def func( self, *args, **kwargs ):
            helix={}
            sheet={}
            fullsheet={}
            sheetnames=[]
            npdb=numpdb.NumPdb( self.pdb_file)
            with open (self.pdb_file) as header:
                headerlines=header.readlines()
                xyz=0
                for lines in headerlines:
                    if lines.startswith  ("HELIX"):
                        ident="%s_%s_%s" % ('HELIX',lines[19:20],lines[8:10])
                        helix[ident]=[int(lines[72:76])]
                        helix[ident].append(lines[8:10])
                        helix[ident].append(lines[31:32])
                        helix[ident].append(int(lines[22:25]))
                        helix[ident].append(lines[34:37])
                    if lines.startswith  ("SHEET"):
                        #print int(lines[34:37])
                        xyz+=1
                        #hurzi=lines.split()
                        #print hurzi[4]
                        ident=lines[12:14],lines[8:10]#"%s_%s_%s.%s" % ("SHEET",lines[21:22],lines[8:10],xyz)
                        sheetlen=(int(lines[34:37])-int(lines[23:26]))+1
                        #print sheetlen
                        sheet[ident]=[sheetlen]#[int(lines[72:76])]
                        sheet[ident].append(lines[8:10])
                        sheet[ident].append(lines[21:22])
                        sheet[ident].append(int(lines[23:26]))
                        sheet[ident].append(int(lines[34:37]))
                        sheet[ident].append(lines[12:14])
                        sheetnames.append(lines[12:14])
            for hurz in sorted(set(sheetnames)):
                go={}
                for key in sorted(sheet.iterkeys()):
                    
                    if key[0]==hurz:
    
                        go[key]=sheet[key]
    
                for i in range (0,len(go)-1):
    
                    if go.values()[i] [2]==go.values()[i+1] [2] and abs(go.values()[i] [4]-go.values()[i+1] [4])<10 :
     
                        sheetlist=[]
                        sheetlist.extend([go.values()[i] [4],go.values()[i] [3],go.values()[i+1] [3],go.values()[i+1] [4]])
                        identfull= "%s_%s_%s" % (hurz,go.values()[i] [2],min(sheetlist))
                        fullsheet[identfull]=[go.values()[i] [2]]
                        fullsheet[identfull].append(min(sheetlist))
                        fullsheet[identfull].append(max(sheetlist))
                        fullsheet[identfull].append(max(sheetlist)-min(sheetlist))
            for key in sorted(fullsheet):
                res1=fullsheet[key][1]
                res2=fullsheet[key][2]
                chain=fullsheet[key][0]
                print key, fullsheet[key]
                test=npdb.copy(resno=[res1,res2],chain=chain)
                seq=test.sequence()
                print seq
                ires1="%s:%s" % (res1-1,chain)
                ires2="%s:%s" % (res2+1,chain)
                diro="%s/%s/%s_%s" % (chain,fullsheet[key][3],res1,res2)
                try:
                    LinkItDensity (self.pdb_file,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=diro)
                except:
                    print 'linikt error'
               # LinkItDensity(self.pdb_file,self.mrc_file,)
class LinkItDensity3( PyTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "mrc_file", type="file", ext="mrc" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "seq", type="str" ),
        _( "resolution", type="float", fixed=True ),
        _( "cutoff", type="float" ),
        _( "max_loops", type="int", range=[0, 500], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it_density.provi"

    def _init( self, *args, **kwargs ):
        if self.res1['resno']>self.res2['resno']:
            self.res1,self.res2=self.res2,self.res1
        self.link_it = LinkIt( 
            self.edited_pdb_file, self.res1, self.res2, self.seq,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("link_it") )
        )
        self.loop_correl = LoopCrosscorrel3 (
            self.mrc_file, self.pdb_file, self.link_it.pdb_linker_file2, self.link_it.txt_file,
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
    def func( self ):
        boxsize=getMrc(self.mrc_file,'nx' )
        originx=abs(getMrc(self.mrc_file,'nxstart' ))
        originy=abs(getMrc(self.mrc_file,'nystart' ))
        originz=abs(getMrc(self.mrc_file,'nzstart' ))
        size=getMrc(self.mrc_file,'xlen' )
        pixelsize=(size/boxsize)
        shx = (originx -(boxsize/2)) * pixelsize
        shy = (originy -(boxsize/2)) * pixelsize
        shz = (originz -(boxsize/2)) * pixelsize
        #print shx,shy,shz
        PdbEdit( 
            self.pdb_file, shift= [shx, shy, shz]
        )    
   
        self.link_it()
        self.loop_correl()
    def _post_exec( self ):
        self._make_correl_json()
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            # pdb_file=self.relpath( self.loop_correl.spider_shift.edited_pdb_file ),
            mrc_file=self.relpath( self.mrc_file ),
            # mrc_file=self.relpath( self.loop_correl.spider_shift.map_shift ),
            cutoff=self.cutoff,
            # box_mrc_file=self.relpath( 
            #     self.loop_correl.spider_reconvert.mrc_file 
            # ),
            # box_ori_mrc_file=self.relpath( 
            #     self.loop_correl.spider_reconvert.mrc_ori_file 
            # ),
            pdb_linker_file3=self.relpath( 
                self.loop_correl.ori_pdb_linker_file3 ),
            linker_correl_file=self.relpath( self.linker_correl_file )
        )
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
        linker_correl_dict = {}
        for k, v in cc_dict.iteritems():
            linker_correl_dict[k] = [ v ] + li_dict[k]
        with open( self.linker_correl_file, "w" ) as fp:
            if compact:
                json.dump( linker_correl_dict, fp, separators=(',',':') )
            else:
                json.dump( linker_correl_dict, fp, indent=4 )

class MultiCutPDB (PyTool):
    args = [
    _( "pdb_file", type="file"),
    _( "mrc_file", type="file"),
    _( "mrc_file2", type="file", ext="mrc" ),
    _( "mrc_file3", type="file", ext="mrc" ),
    _( "mrc_file4", type="file", ext="mrc" ),
    _( "mrc_file5", type="file", ext="mrc" ),
    _( "mrc_file6", type="file", ext="mrc" ),
    _( "mrc_file7", type="file", ext="mrc" ),
    _( "resolution", type="float", range=[1, 20], step=0.1 ),
    _( "resolution2", type="float", fixed=True ),
    _( "resolution3", type="float", fixed=True ),
    _( "resolution4", type="float", fixed=True ),
    _( "resolution5", type="float", fixed=True ),
    _( "resolution6", type="float", fixed=True ),
    _( "resolution7", type="float", fixed=True ),
    _( "cutoff", type="float" )
        ]
    def func( self, *args, **kwargs ):
        
        npdb=numpdb.NumPdb( self.pdb_file)

        for z, numa in enumerate (npdb.iter_chain()):
            cha=numa.get('chain')[0]
            rel=numa.get('resno')[-1]
            print rel
            print cha
            if cha in ('P','h','i'):
                os.mkdir(cha)
                for x in range (12,20,1):
                    ch=numa.get('chain')[0]
                    gum= "%s/%s" % (cha,x)
                    os.mkdir(gum)
                    for i,numa in enumerate (npdb.iter_resno(chain=cha)):
                       
                        
                       
                        res1=numa.get('resno')[0]
                        if res1 % 7 == 0 :#-->jeder 7. loop
                            res2=res1+x-1
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
                               # start = timeit.timeit()
                                #print 'richtige?',numa.sequence()
                                try:
                                    MultiLinkItDensity (name2,self.mrc_file,self.mrc_file2,self.mrc_file3,self.mrc_file4,self.mrc_file5,self.mrc_file6,self.mrc_file7,ires1,ires2,seq,self.resolution,self.resolution2,self.resolution3,self.resolution4,self.resolution5, self.resolution6,self.resolution7,self.cutoff,output_dir=di)
                                except:
                                    print 'linikt error'
                                #    neuen loop suchen
                #break
class MultiLinkItDensity( PyTool, ProviMixin ):
    """Please choose a cutoff  and denote the resolution of your map. Then define the stem-residues (e.g.: 10:A, 16:A) and provide the missing sequence in 1-letter code (ACDEF). """
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "mrc_file", type="file", ext="mrc" ),
        _( "mrc_file2", type="file", ext="mrc" ),
        _( "mrc_file3", type="file", ext="mrc" ),
        _( "mrc_file4", type="file", ext="mrc" ),
        _( "mrc_file5", type="file", ext="mrc" ),
        _( "mrc_file6", type="file", ext="mrc" ),
        _( "mrc_file7", type="file", ext="mrc" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "seq", type="str" ),
        _( "resolution", type="float", fixed=True ),
        _( "resolution2", type="float", fixed=True ),
        _( "resolution3", type="float", fixed=True ),
        _( "resolution4", type="float", fixed=True ),
        _( "resolution5", type="float", fixed=True ),
        _( "resolution6", type="float", fixed=True ),
        _( "resolution7", type="float", fixed=True ),
        _( "cutoff", type="float" ),
        _( "max_loops", type="int", range=[0, 500], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it_density.provi"
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
        self.loop_correl2 = LoopCrosscorrel(
            self.mrc_file2, self.pdb_file,
            self.link_it.pdb_linker_file2,
            self.link_it.txt_file,
            self.res1, self.res2, len(self.seq),
            self.resolution2,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("loop_correl2 "),
                max_loops=self.max_loops,
            )
        )
        self.loop_correl3 = LoopCrosscorrel(
            self.mrc_file3, self.pdb_file,
            self.link_it.pdb_linker_file2,
            self.link_it.txt_file,
            self.res1, self.res2, len(self.seq),
            self.resolution3,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("loop_correl3 "),
                max_loops=self.max_loops,
            )
        )
        self.loop_correl4 = LoopCrosscorrel(
            self.mrc_file4, self.pdb_file,
            self.link_it.pdb_linker_file2,
            self.link_it.txt_file,
            self.res1, self.res2, len(self.seq),
            self.resolution4,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("loop_correl4 "),
                max_loops=self.max_loops,
            )
        )
        self.loop_correl5 = LoopCrosscorrel(
            self.mrc_file5, self.pdb_file,
            self.link_it.pdb_linker_file2,
            self.link_it.txt_file,
            self.res1, self.res2, len(self.seq),
            self.resolution5,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("loop_correl5 "),
                max_loops=self.max_loops,
            )
        )
        self.loop_correl6 = LoopCrosscorrel(
            self.mrc_file6, self.pdb_file,
            self.link_it.pdb_linker_file2,
            self.link_it.txt_file,
            self.res1, self.res2, len(self.seq),
            self.resolution6,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("loop_correl6 "),
                max_loops=self.max_loops,
            )
        )
        self.loop_correl7 = LoopCrosscorrel(
            self.mrc_file7, self.pdb_file,
            self.link_it.pdb_linker_file2,
            self.link_it.txt_file,
            self.res1, self.res2, len(self.seq),
            self.resolution7,
            **copy_dict(
                kwargs, run=False, output_dir=self.subdir("loop_correl7 "),
                max_loops=self.max_loops,
            )
        )
        self.output_files += list( itertools.chain(
            self.link_it.output_files,
            self.loop_correl.output_files,
        ))

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
        print originx
        if originx!=0 and xorg==0 :
            shx = (originx -(boxsizex/2)) * pixelsize
            shy = (originy -(boxsizey/2)) * pixelsize
            shz = (originz -(boxsizez/2)) * pixelsize
        else:
            shx = ((xorg  / pixelsize)-(boxsizex/2))*pixelsize
            shy = ((yorg/ pixelsize)-(boxsizey/2))*pixelsize
            shz = ((zorg / pixelsize)-(boxsizez/2))*pixelsize
        PdbEdit(
            self.pdb_file, shift=[0,0,0]#[ shx, shy, shz ]
        )

        self.link_it()
        self.loop_correl()
        self.loop_correl2()
        self.loop_correl3()
        self.loop_correl4()
        self.loop_correl5()
        self.loop_correl6()
        self.loop_correl7()
        print len(self.seq)
    def _post_exec( self ):
        self._make_correl_json()
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            # pdb_file=self.relpath( self.loop_correl.spider_shift.edited_pdb_file ),
            mrc_file=self.relpath( self.mrc_file ),
            # mrc_file=self.relpath( self.loop_correl.spider_shift.map_shift ),
            cutoff=self.cutoff,
            # box_mrc_file=self.relpath(
            #     self.loop_correl.spider_reconvert.mrc_file
            # ),
            # box_ori_mrc_file=self.relpath(
            #     self.loop_correl.spider_reconvert.mrc_ori_file
            # ),
            pdb_linker_file3=self.relpath(
                self.loop_correl.ori_pdb_linker_file3 ),
            linker_correl_file=self.relpath( self.linker_correl_file )
        )
        
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
        linker_correl_dict = {}
        for k, v in cc_dict.iteritems():
            linker_correl_dict[k] = [ v ] + li_dict[k]
        with open( self.linker_correl_file, "w" ) as fp:
            if compact:
                json.dump( linker_correl_dict, fp, separators=(',', ':') )
            else:
                json.dump( linker_correl_dict, fp, indent=4 )
    
    #end = timeit.timeit()
    #print "zEdIT", end - start
class LinkItList (PyTool):
    args = [
        _( "linklist", type="file", ext="txt" ),
        _( "dataset_dir", type="dir" )
        ]
    def func( self ):
        with open(self.linklist,'r') as linki,open('fails.txt','w') as fail:
            gul=linki.readlines()
            for line in gul[1:-1]:
                pdb=line.split( )[0]
                as1=line.split( )[1]
                as2=line.split( )[2]
                seq=line.split( )[3]
                res1='%s:%s' % (as1,'A')
                res2='%s:%s' % (as2,'A')
                pdbfile='%s/%s.%s' % (self.dataset_dir,pdb,'pdb')
                outdir='%s_%s_%s' % (pdb,as1,as2)
                print pdb,as1,as2, seq[1:-1]
                try:
                    LinkIt(pdbfile,res1,res2,seq[1:-1],maxloops=500,output_dir=self.subdir(outdir))
                except:
                    print 'linkiterror'
                    fail.write(outdir)
                    fail.write('\n')
class PdbSequence (PyTool):
    args = [
        _( "dataset_dir", type="dir" )
         ]
    def func( self ):
        hurz=[h for h in os.listdir(self.dataset_dir) if any ([h.endswith('.pdb')]) ]
        for i in hurz:
            name="%s_%s.%s" % (i,'sequence','txt')
            pdbname=os.path.join(self.dataset_dir,i)
            npdb=numpdb.NumPdb(pdbname)
            out=os.path.join('seq_sali',name)
            with open (out,'w') as seqfile:
                for z, numa in enumerate (npdb.iter_chain()):
                    cha=numa.get('chain')[0]
                    seq=numa.sequence()
                    seqclean=seq.replace("?","")
                    seqfile.write('>')
                    seqfile.write(cha)
                    seqfile.write('\n')
                    seqfile.write(seqclean)
                    seqfile.write('\n')
            seqfile.close()
            