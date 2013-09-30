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


import provi_prep as provi
from spider import LoopCrosscorrel
from spider import MrcHeaderPrinter, MrcHeader, mrc_header
from pdb import PdbEdit, CutpdbSSE 

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "linker" ) 


LINKIT_DIR = os.environ.get("LINKIT_DIR", "")
LINKIT_CMD = os.path.join( LINKIT_DIR, "Link_It_dos2n.exe" )



class LinkerTest( PyTool ):
    args = [
        _( "linker_txt", type="file", ext="txt" )
    ]
    out = [
        _( "linker_json", file="{linker_txt.stem}.json" )
    ]
    def func( self ):
        self._make_linker_json( compact=True )
    def _make_linker_json( compact=False ):
        linker_dict = {}
        with open( self.linker_txt, "r" ) as fp:
            x = fp.next()
            n = fp.next()
            for i, d in enumerate( iter_stride( fp, 3 ), start=1 ):
                linker_dict[ i ] = [ float(d[0]), float(d[1]), d[2].strip() ]
        with open( self.linker_json, "w" ) as fp:
            if compact:
                json.dump( linker_dict, fp, separators=(',',':') )
            else:
                json.dump( linker_dict, fp, indent=4 )


class LinkIt( CmdTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "seq", type="text" )
    ]
    out = [
        _( "bin_file", file="{pdb_file.stem}_linker.bin" ),
        _( "txt_file", file="{pdb_file.stem}_linker.txt" ),
        _( "pdb_linker_file", file="{pdb_file.stem}_linker.pdb" ),
        _( "pdb_linker_file2", file="{pdb_file.stem}_linker2.pdb" ),
        _( "pdb_linker_file3", file="{pdb_file.stem}_linker3.pdb" ),
        _( "kos_file", file="{pdb_file.stem}_kos.txt" ),
        _( "json_file", file="{pdb_file.stem}_linker.json" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it.provi"
    def _init( self, *args, **kwargs ):
        self.cmd = [ "wine", LINKIT_CMD, self.kos_file, self.bin_file, "t" ]
    def _pre_exec( self ):
        self._make_kos_file()
    def _post_exec( self ):
        self._fix_linker_pdb( self.pdb_linker_file2 )
        self._fix_linker_pdb( self.pdb_linker_file3, atoms_only=True )
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
                coords = npdb.get( 'xyz', **sele )
                fp.write( "%s\n" % "\n".join(map( str, coords[0] ) ) )
            fp.write( "%s\n" % self.seq )
            for sele in ( self.res1, self.res2 ):
                fp.write( "%s\t%s\n" % ( sele.get("chain", " "), sele["resno"] ) )
    def _fix_linker_pdb( self, output_file, atoms_only=False, stems=True ):
        backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
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
                        atom_i += 1
                        fp_out.write( line )
                        continue
                    if not atoms_only:
                        fp_out.write( line )
    def _make_linker_json( self, compact=False ):
        linker_dict = {}
        with open( self.txt_file, "r" ) as fp:
            x = fp.next()
            n = fp.next()
            for i, d in enumerate( iter_stride( fp, 4 ), start=1 ):
                linker_dict[ i ] = [ float(d[0]), float(d[1]), str(d[2].strip()),str(d[3].strip()) ]
                #print d
        with open( self.json_file, "w" ) as fp:
            if compact:
                json.dump( linker_dict, fp, separators=(',',':') )
            else:
                json.dump( linker_dict, fp, indent=4 )

def getMrc( mrc_file, param ):
    header = mrc_header( mrc_file )
    for name, value in zip(header._fields, header):
        if name==param:
            return value    

class LinkItDensity( PyTool, ProviMixin ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "mrc_file", type="file", ext="mrc" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "seq", type="text" ),
        #_( "pixelsize", type="slider", range=[1, 10], fixed=True ),
        _( "resolution", type="slider", range=[1, 10], fixed=True ),
        #_( "boxsize", type="slider", range=[1, 500], fixed=True ),
        #_( "originx", type ="slider", range=[-500,500]),
        #_( "originy", type ="slider", range=[-500,500]),
        #_( "originz", type ="slider", range=[-500,500]),
        _( "cutoff", type="float", default=5 ),
        _( "max_loops", type="slider", range=[0, 200], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it_density.provi"

    def _init( self, *args, **kwargs ):
        self.link_it = LinkIt( 
            self.edited_pdb_file, self.res1, self.res2, self.seq,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("link_it") )
        )
        self.loop_correl = LoopCrosscorrel(
            self.mrc_file, self.pdb_file, self.link_it.pdb_linker_file2, 
            self.res1, self.res2, len(self.seq),
            self.resolution,            
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("loop_correl"),
                max_loops=self.max_loops,
            )
        )
        self.output_files = list( itertools.chain(
            self.link_it.output_files, 
            self.loop_correl.output_files,
        ))
    def func( self ):
        boxsize=getMrc(self.mrc_file,'nx' )
        originx=abs(getMrc(self.mrc_file,'xorg' ))
        originy=abs(getMrc(self.mrc_file,'yorg' ))
        originz=abs(getMrc(self.mrc_file,'zorg' ))
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
            mrc_file=self.relpath( self.mrc_file ),
            cutoff=self.cutoff,
            box_mrc_file=self.relpath( 
                self.loop_correl.spider_reconvert.mrc_file 
            ),
            box_ori_mrc_file=self.relpath( 
                self.loop_correl.spider_reconvert.mrc_ori_file 
            ),
            pdb_linker_file3=self.relpath( self.link_it.pdb_linker_file3 ),
            linker_correl_file=self.relpath( self.linker_correl_file )
        )
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



class LinkItDensityjason( PyTool, ProviMixin ):
    args = [
        _( "json_file", type="file", ext=".json" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "mrc_file", type="file", ext="mrc" ),
        _( "xml_file", type="file", ext="xml" ),
        #_( "res1", type="sele" ),
        #_( "res2", type="sele" ),
        #_( "seq", type="text" ),
        #_( "pixelsize", type="slider", range=[1, 10], fixed=True ),
        #_( "resolution", type="slider", range=[1, 10], fixed=True ),
        #_( "boxsize", type="slider", range=[1, 500], fixed=True ),
        #_( "originx", type ="slider", range=[-500,500]),
        #_( "originy", type ="slider", range=[-500,500]),
        #_( "originz", type ="slider", range=[-500,500]),
        _( "cutoff", type="float", default=5 ),
        _( "max_loops", type="slider", range=[0, 200], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it_density.provi"

    def _init( self, *args, **kwargs ):
        with open(self.json_file,"r") as hz:
            js=json.load(hz)
            print js["3"]
            res1=js["1"]
            res2=js["2"]
            seq=js["3"]
            chain=js["4"]
            stem1=str(res1)+":"+chain
            stem2=str(res2)+":"+chain
            print stem1
        with open (self.xml_file, "r") as rs:
            tree = ET.parse(rs)
            root = tree.getroot()
            #for child in root:
            #    print child.tag, child.attrib
            resu=float(root.find("processing/reconstruction/resolutionByAuthor").text)
            #resu=doc.find("resolutionByAuthor")
            print resu
        self.link_it = LinkIt( 
            self.edited_pdb_file,stem1,stem2, seq,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("link_it") )
        )
        self.loop_correl = LoopCrosscorrel(
            self.mrc_file, self.pdb_file, self.link_it.pdb_linker_file2, 
            stem1, stem2, len(seq),
            resu,            
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("loop_correl"),
                max_loops=self.max_loops,
            )
        )
        self.output_files = list( itertools.chain(
            self.link_it.output_files, 
            self.loop_correl.output_files,
        ))
    def func( self ):
        boxsize=getMrc(self.mrc_file,'nx' )
        originx=abs(getMrc(self.mrc_file,'nxstart'))#'xorg' ))
        originy=abs(getMrc(self.mrc_file,'nystart'))#'yorg' ))
        originz=abs(getMrc(self.mrc_file,'nzstart'))#'zorg' ))
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
            mrc_file=self.relpath( self.mrc_file ),
            cutoff=self.cutoff,
            box_mrc_file=self.relpath( 
                self.loop_correl.spider_reconvert.mrc_file 
            ),
            box_ori_mrc_file=self.relpath( 
                self.loop_correl.spider_reconvert.mrc_ori_file 
            ),
            pdb_linker_file3=self.relpath( self.link_it.pdb_linker_file3 ),
            linker_correl_file=self.relpath( self.linker_correl_file )
        )
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
        _( "max_loops", type="slider", range=[0, 200], default=100 )
    ]
    out = [
        _( "linker_correl_file", file="linker_correl.json" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    def func( self ):
        
        for fn in os.listdir(self.dataset_dir):
            print fn
            for files in os.listdir(self.dataset_dir+"/"+fn+"/pieces"):
                if files.endswith(".json"):
                    jason=self.dataset_dir+"/"+fn+"/pieces/"+files
                    pdb=self.dataset_dir+"/"+fn+"/"+files [0:4]+".pdb"
                    emmap=self.dataset_dir+"/"+fn+"/map/"+fn[0:3].lower()+"_"+fn[4:8]+".map"
                    xml_file=self.dataset_dir+"/"+fn+"/header/"+fn.lower()+".xml"
                    #print pdb
                    print emmap
                    print xml_file
                    LinkItDensityjason (
                        jason, pdb, emmap, xml_file,  output_dir=self.dataset_dir+"/"+fn+"/pieces/"+files[:-4]                   
                    )
                    #output verzeichnis init functio aufteilung?
                    
