from __future__ import with_statement
from __future__ import division

import os
import json
import itertools
from string import Template

from utils import copy_dict
from utils.tool import _, _dir_init, PyTool, CmdTool, ScriptMixin, ProviMixin
from utils.numpdb import NumPdb, numsele

from pdb import PdbSplit
from pdb import PdbEdit

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "spider" )
SPIDER_CMD = "spider" 


# 2010 Cryo-EM Modeling Challenge: http://ncmi.bcm.edu/challenge


class Spider( CmdTool, ScriptMixin ):
    args = [
        _( "script_file", type="file", ext="spi" )
    ]
    script_tmpl = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, *args, **kwargs ):
        if script_file=="__tmpl__":
            self.script_file = self.outpath( self.script_tmpl )
        script_ext = "spi"
        data_ext = "cpv"
        # spider spi/cpv @box
        self.cmd = [
            SPIDER_CMD, 
            "%s/%s" % ( script_ext, data_ext ), 
            "@%s" % self.relpath( self.script_file, no_ext=True )
        ]

class SpiderShift( Spider ):
    """Shifts and converts mrc file and pdb"""
    args = [
        _( "mrc_file", type="file", ext="map" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "pixelsize", type="slider", range=[1, 10], fixed=True ),
        _( "boxsize", type="slider", range=[1, 500], fixed=True ),
        _( "originx", type ="slider", range=[-500,500]),
        _( "originy", type ="slider", range=[-500,500]),
        _( "originz", type ="slider", range=[-500,500]),
        #_( "shift", type="float", nargs=3, default=None,
        #   metavar=("X", "Y", "Z") ),
    ]
    out = [
        _( "map_shift", file="shift.mrc" ),
        _( "edited_pdb_file", file="edited.pdb" )
    ]  
    script_tmpl = "shift.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderShift, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file(    
            mrc_file=self.relpath( self.mrc_file ),
            pixelsize=self.pixelsize,
            boxsize=self.boxsize,
            originx=self.originx,
            originy=self.originy,
            originz=self.originz,
        )
        
        shx = (self.originx -(self.boxsize/2)) * self.pixelsize
        shy = (self.originy -(self.boxsize/2)) * self.pixelsize
        shz = (self.originz -(self.boxsize/2)) * self.pixelsize
        PdbEdit( 
            self.pdb_file, shift= [shx, shy, shz]
        )    

class SpiderConvert( Spider ):
    """Simple tool that converts mrc files to the spider format"""
    args = [
        _( "mrc_file", type="file", ext="mrc" )
    ]
    out = [
        _( "map_file", file="mapupload.cpv" )
    ]
    script_tmpl = "convert.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderConvert, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file( 
            mrc_file=self.relpath( self.mrc_file )
        )

class SpiderPdbBox( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "result_file", type="file" ),
        _( "boxsize", type="slider", range=[1, 500], fixed=True )
    ]
    out = [
        _( "edited_pdb_file", file="edited.pdb" )
    ]  
    def func( self ):
        with open( self.result_file, "r") as rf:
            content = rf.readlines ()
            al = content[1].strip ()
            rs = al.split ()
            bs = float (rs [2])
            ol1 = float (rs [3])
            ol2 = float (rs [4])
            ol3 = float (rs [5])
            ps = float (rs [7])
            x =  (ol1 - (self.boxsize/2)) * ps - 1 
            y =  (ol2 - (self.boxsize/2)) * ps - 1
            z =  (ol3 - (self.boxsize/2)) * ps - 1
            #print [x,y,z, bs]
        PdbEdit (self.pdb_file, box = [ x, y, z, bs, bs, bs ] )
    
class SpiderDeleteFilledDensities( Spider ):
    args = [
        _( "map_file", type="file", ext="cpv" , help= "boxil.cpv" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "result_file", type="file" , ext="cpv", help="ergebnisse.cpv"),
        _( "pixelsize", type="slider", range=[1, 10], fixed=True ),
        _( "resolution", type="slider", range=[1, 10], 
            fixed=True, help="of the map_file" ),
        _( "boxsize", type="slider", range=[1, 500], fixed=True , help = "of the input map")
    ]
    out = [
        _( "empty_map_file", file="usermap.cpv" )
    ]
    script_tmpl = "delete_filled_densities.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderDeleteFilledDensities, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            pdb_file=self.relpath( self.pdb_file ),
            result_file=self.relpath( self.result_file, no_ext=True ),
            pixelsize=self.pixelsize,
            resolution=self.resolution,
            boxsize=self.boxsize,
            tmp_dir=self.relpath( self.output_dir ) + os.sep
        )



class SpiderBox( Spider ):
    args = [
        _( "map_file", type="file", ext="cpv" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "res1", type="sele", help="resno:chain, i.e. 10:A" ),
        _( "res2", type="sele" ),
        _( "length", type="slider", range=[1, 30] ),
        _( "pixelsize", type="slider", range=[1, 10], fixed=True ),
        _( "resolution", type="slider", range=[1, 10], 
            fixed=True, help="of the map_file" )
    ]
    out = [
        _( "var_file", file="variables.cpv" ),
        _( "box_file", file="ergebnisse.cpv" ),
        _( "box_map_file", file="boxil.cpv" )
    ]
    script_tmpl = "box.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderBox, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        coords1, coords2 = self._get_coords( 
            self.pdb_file, self.res1, self.res2 
        )
        self._make_variables_file(
            coords1, coords2, self.length, self.pixelsize, self.resolution
        )
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            var_name=self.relpath( self.var_file, no_ext=True )
        )
    def _get_coords( self, pdb_file, res1, res2 ):
        npdb = NumPdb( pdb_file )
        return npdb.center( **res1 ), npdb.center( **res2 )
    def _make_variables_file( self, coords1, coords2, length, 
                              pixelsize, resolution ):
        variables = "1 9 %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %i %4.2f %4.2f" % (
            coords1[0], coords1[1], coords1[2],
            coords2[0], coords2[1], coords2[2],
            length, pixelsize, resolution
        )
        with open( self.var_file, "w" ) as fp:
            fp.write( variables )


class SpiderReConvert( Spider ):
    args = [
        _( "box_file", type="file", ext="cpv" ),
        _( "map_file", type="file", ext="cpv" ),
        _( "box_map_file", type="file", ext="cpv" ),
        _( "ori_map_file", type="file", ext="cpv" )
    ]
    out = [
        _( "mrc_file", file="reconvert.mrc" ),
        _( "mrc_ori_file", file="reconvertori.mrc" )
    ]
    script_tmpl = "recon.spi"
    def _init( self, *args, **kwargs ):
        super(SpiderReConvert, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file( 
            box_name=self.relpath( self.box_file, no_ext=True ),
            map_name=self.relpath( self.map_file, no_ext=True ),
            box_map_name=self.relpath( self.box_map_file, no_ext=True ),
            ori_map_name=self.relpath( self.ori_map_file, no_ext=True )
        )


class SpiderCrosscorrelation( Spider ):
    args = [
        _( "map_file", type="file", ext="cpv" ),
        _( "box_map_file", type="file", ext="cpv" ),
        _( "box_file", type="file", ext="cpv" ),
        _( "loop_file", type="file", ext="pdb" ),
        _( "max_loops", type="slider", range=[0, 200], default=100 )
    ]
    out = [
        _( "loop_dir", dir="loops" ),
        _( "crosscorrel_file", file="crosscorrelation.cpv" ),
        _( "crosscorrel_json", file="crosscorrelation.json" )
    ]
    script_tmpl = "crosscorrelation.spi"
    def _init( self, *args, **kwargs ):
        if not self.max_loops:
            self.max_loops = False
        super(SpiderCrosscorrelation, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._split_loop_file()
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            box_map_name=self.relpath( self.box_map_file, no_ext=True ), 
            box_name=self.relpath( self.box_file, no_ext=True ), 
            loop_dir=self.relpath( self.loop_dir ) + os.sep,
            max_loops=self.max_loops or 999
        )
    def _post_exec( self ):
        self._make_crosscorrel_json( compact=True )
    def _split_loop_file( self ):
        PdbSplit( 
            self.loop_file, output_dir=self.loop_dir, backbone_only=True, 
            max_models=self.max_loops, resno_ignore=[ 1000, 2000 ], zfill=3
        )
    def _make_crosscorrel_json( self, compact=False ):
        crosscorrel_dict = {}
        with open( self.crosscorrel_file, "r" ) as fp:
            for line in fp:
                d = line.split()
                if len(d)==3:
                    crosscorrel_dict[ int(d[0]) ] = float( d[2] )
        with open( self.crosscorrel_json, "w" ) as fp:
            if compact:
                json.dump( crosscorrel_dict, fp, separators=(',',':') )
            else:
                json.dump( crosscorrel_dict, fp, indent=4 )


class LoopCrosscorrel( PyTool ):
    args = [
        _( "mrc_file", type="file", ext="mrc" ),
        _( "pdb_file", type="file", ext="pdb" ),
        _( "loop_file", type="file", ext="pdb" ),
        _( "res1", type="sele" ),
        _( "res2", type="sele" ),
        _( "length", type="slider", range=[1, 30] ),
        _( "pixelsize", type="slider", range=[1, 10], fixed=True ),
        _( "resolution", type="slider", range=[1, 10], fixed=True ),
        _( "boxsize", type="slider", range=[1, 500], fixed=True ),
        _( "originx", type ="slider", range=[-500,500]),
        _( "originy", type ="slider", range=[-500,500]),
        _( "originz", type ="slider", range=[-500,500]),
        _( "max_loops", type="slider", range=[0, 200], default=100 )
    ]
    out = [
        _( "cropped_pdb", file="cropped.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        self.spider_shift = SpiderShift(
            self.mrc_file,self.pdb_file,self.pixelsize, self.boxsize,
            self.originx,self.originy,self.originz,
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("shift") 
            )
        )
        self.spider_convert = SpiderConvert( 
            self.spider_shift.map_shift, 
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("convert") 
            )
        )
        self.spider_box = SpiderBox(
            self.spider_convert.map_file, 
            self.spider_shift.edited_pdb_file, self.res1, self.res2, 
            self.length, self.pixelsize, self.resolution,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("box") )
        )
        self.pdb_box =  SpiderPdbBox(
           self.spider_shift.edited_pdb_file,
           self.spider_box.box_file,
           self.boxsize,
           **copy_dict( kwargs, run=False, output_dir=self.subdir("pdbbox"))
        )
        self.spider_delete_filled_densities = SpiderDeleteFilledDensities( 
            self.spider_box.box_map_file, self.pdb_box.edited_pdb_file, self.spider_box.box_file,
            self.pixelsize, self.resolution,self.boxsize,
            **copy_dict( 
                kwargs, run=False, 
                output_dir=self.subdir("delete_filled_densities") 
            )
        )   


        self.spider_reconvert = SpiderReConvert(
            self.spider_box.box_file,
            self.spider_convert.map_file,
            self.spider_box.box_map_file,
            self.spider_convert.map_file,
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("reconvert") 
            )
        )
        self.spider_crosscorrelation = SpiderCrosscorrelation(
            self.spider_convert.map_file, 
            self.spider_box.box_map_file, 
            self.spider_box.box_file, 
            self.loop_file,
            **copy_dict( 
                kwargs, run=False, output_dir=self.subdir("crosscorrelation"),
                max_loops=self.max_loops
            )
        )
        self.output_files.extend( list( itertools.chain(
            self.spider_convert.output_files, 
            self.spider_delete_filled_densities.output_files,
            self.spider_box.output_files,
            self.spider_reconvert.output_files,
            self.spider_crosscorrelation.output_files
        )))
    def func( self ):
        self._crop_pdb()
        self.spider_shift()
        self.spider_convert()
        self.spider_box()
        self.spider_pdb_box ()
        self.spider_delete_filled_densities()
        self.spider_reconvert()
        self.spider_crosscorrelation()
    def _crop_pdb( self ):
        npdb = NumPdb( self.pdb_file )
        npdb.write( 
            self.cropped_pdb, 
            chain=self.res1["chain"], 
            resno=[ self.res1["resno"]+1, self.res2["resno"]-1 ],
            invert=True
        )



