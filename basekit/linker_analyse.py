from __future__ import with_statement
from __future__ import division
                                                                                                                                                                                                                                                                                                                                                                                                                        

import os
import itertools
import operator
import json
import collections
import xml.etree.ElementTree as ET

import utils.path
from utils import copy_dict, iter_stride
from utils.tool import _, _dir_init, CmdTool, PyTool, ProviMixin
from utils.numpdb import NumPdb, numsele
from utils.mrc import get_mrc_header, getMrc
import numpy as np
import time
import provi_prep as provi
from spider_analyse import LoopCrosscorrel,LoopCrosscorrel2,LoopCrosscorrel3
from pdb import PdbEdit, SplitPdbSSE, LoopDelete ,PdbSplit,get_tree
import utils.numpdb as numpdb

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "linker" ) 

def LINKIT_DIR():
    return os.environ.get("LINKIT_DIR", "")

def LINKIT_CMD():
    return os.path.join( LINKIT_DIR(), "Link_It_dos2n.exe" )



class LinkIt( CmdTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "seq", type="str" )
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
    #t1=time.time()
    #print t1
    def _init( self, *args, **kwargs ):
        if self.res1['resno']>self.res2['resno']:
            self.res1, self.res2 = self.res2, self.res1
        print 'dtart linkit',time.time()    
        self.cmd = [ "wine", LINKIT_CMD(), self.kos_file, self.bin_file, "t" ]
        
        
    def _pre_exec( self ):
        self._make_kos_file()
    def _post_exec( self ):
        
        print 'end linkit', time.time()
        self._fix_linker_pdb( self.pdb_linker_file2 )
        print '1.fix', time.time()
        self._fix_linker_pdb( self.pdb_linker_file3, atoms_only=True )
        print '2.fix', time.time()
        self._split_loop_file()
        print '2.fix', time.time()
        self._make_linker_json( compact=True )
        print 'make linker_json', time.time()
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            pdb_linker_file3=self.relpath( self.pdb_linker_file3 ),
            json_file=self.relpath( self.json_file )
        
        )
        print 'make provi',time.time()
        #print t3-t1
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
                coords = npdb.get( 'xyz', **sele )
                fp.write( "%s\n" % "\n".join(map( str, coords[0] ) ) )
            fp.write( "%s\n" % self.seq )
            for sele in ( self.res1, self.res2 ):
                fp.write( 
                    "%s %s\n" % ( sele.get("chain") or " ", sele["resno"] )
                )
    def _fix_linker_pdb( self, output_file, atoms_only=False, stems=True ):
        backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
        chain = self.res1['chain']
        print chain
        #print time.time()
        with open( self.pdb_linker_file, "r" ) as fp:
            with open( output_file, "w" ) as fp_out:
                for i, line in enumerate( fp ):
                    if line.startswith("MODEL"):
                        atom_i = 1
                    if line.startswith("ATOM"):
                        line = line = line[0:6] + ( "% 5i" % atom_i ) + line[11:]
                        if line[22]=="X":
                            if not stems:
                                continue
                            if line[12:16] not in backbone:
                                continue
                            tag = "1000 " if line[24]==" " else "2000 "
                            line = line[0:17] + "GLY" + line[20:22] + tag + line[27:]
                        else:
                            resnew= int(line[22:26])+int(self.res1['resno'])
        
                            resnewp="%4s" % resnew
                            line = line = line[0:21]+chain + resnewp + line[26:]
                        atom_i += 1
                        fp_out.write( line )
                        continue
                    if not atoms_only:
                        fp_out.write( line )
        #print time.time()
    def _split_loop_file( self ):
        PdbSplit( 
            self.pdb_linker_file2, output_dir=self.loop_dir, backbone_only=True, 
             resno_ignore=[ 1000, 2000 ], zfill=3
        )                    
    def _find_clashes ( self ):
        backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
        npdb=NumPdb( self.pdb_file ,{"backbone_only": True})
        protein=get_tree(npdb['xyz'])
       
        clashing_models=[]
        
        for  fn in os.listdir(self.loop_dir):
    
            if fn.endswith(".pdb"):
                lf=os.path.join(self.loop_dir,fn)
                npdb2=NumPdb( lf,
            {"backbone_only": True} )
                loops=get_tree(npdb2['xyz'])
                k=loops.query_ball_tree(protein, 3)
                g = [x for x in k if x != []]
    
                f=itertools.chain(*g)
                clashatoms=sorted(set(list(f)))
                clashes=[]
           
                for i in clashatoms:
                    e=npdb.get('resno')[i]
                    if e not in (self.res1['resno'],self.res2['resno']):
                        clashes.append(e)
                if len(clashes)!=0:
                    model=fn.split('_')[0]
                    clashing_models.append(float(model))
    
        return clashing_models
    def _make_linker_json( self, compact=False ):
        linker_dict = {}
        #clashing_models =self._find_clashes()
        #
        #with open( self.txt_file, "r" ) as fp:
        #    fp.next()
        #    fp.next()
        #    for i, d in enumerate( iter_stride( fp, 4 ), start=1 ):
        #        if i in clashing_models :
        #            flag=1
        #        else:
        #            flag=0
        #        linker_dict[ i ] = [ float(d[0]), float(d[1]), str(d[2].strip()),str(d[3].strip()),flag ]

        with open( self.json_file, "w" ) as fp:
            if compact:
                json.dump( linker_dict, fp, separators=(',',':') )
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
        self.loop_correl = LoopCrosscorrel(
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
            for x in range (20,35,1):
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
                            #print 'richtige?',numa.sequence()
                            try:
                                LinkItDensity (name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=di)
                            except:
                                print 'linikt error'
                            #    neuen loop suchen
            break
                            

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
            cha=numa.get('chain')[0]
            rel=numa.get('resno')[-1]
            #if cha!='0':
            #if cha in ('1') : 
            print rel
            print cha
            try:
                os.mkdir(cha)
            except:
                continue
            for x in range (30,36,1):
                ch=numa.get('chain')[0]
                gum= "%s/%s" % (cha,x)
                os.mkdir(gum)
                for i,numa in enumerate (npdb.iter_resno(chain=cha)):
                    #numa.next()
                    #numa.next()
                    #numa.next()
                    res1=numa.get('resno')[0]
                    if res1 % 3 == 0 :#-->jeder 7. loop
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
                                    LinkItDensity2 (name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=di)
                                except:
                                    print 'linikt error'
                                #    neuen loop suchen
                            except:
                                continue
                    #else:
                    #    continue
                    #        
class LinkItDensity2( PyTool, ProviMixin ):
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
        self.loop_correl = LoopCrosscorrel2 (
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
            for x in range (34,36,1):
                ch=numa.get('chain')[0]
                gum= "%s/%s" % (cha,x)
                
                os.mkdir(gum)
                if ch=='A':
                    
                    for i,numa in enumerate (npdb.iter_resno(chain=cha)):
                        #numa.next()
                        #numa.next()
                        #numa.next()
                        res1=numa.get('resno')[0]
                        if res1 % 3 == 0 :#-->jeder 7. loop
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
                                        LinkItDensity3 (name2,self.mrc_file,ires1,ires2,seq,self.resolution, self.cutoff,output_dir=di)
                                    except:
                                        print 'linikt error'
                                except:
                                    continue
                                #    neuen loop suchen
                        
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

