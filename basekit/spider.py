from __future__ import with_statement
from __future__ import division

import os
import json
import itertools
from string import Template

from utils import copy_dict
from utils.tool import PyTool, CmdTool, ScriptMixin, ProviMixin
from utils.numpdb import NumPdb, numsele

from pdb import PdbSplit


DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "spider" )
SPIDER_CMD = "spider" 





class Spider( CmdTool, ScriptMixin ):
    args = [
        { "name": "script_file", "type": "file", "ext": "spi" }
    ]
    script_tmpl = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, script_ext="spi", data_ext="cpv", **kwargs ):
        if script_file=="__tmpl__":
            script_file = self.outpath( self.script_tmpl )
        self.script_file = self.abspath( script_file )
        # spider spi/cpv @box
        self.cmd = [
            SPIDER_CMD, 
            "%s/%s" % ( script_ext, data_ext ), 
            "@%s" % self.relpath( self.script_file, no_ext=True )
        ]



class SpiderConvert( Spider ):
    """Simple tool that converts mrc files to the spider format"""
    args = [
        { "name": "mrc_file", "type": "file", "ext": "mrc" }
    ]
    script_tmpl = "convert.spi"
    def _init( self, mrc_file, **kwargs ):
        self.mrc_file = self.abspath( mrc_file )
        self.map_file = self.outpath( "mapupload.cpv" )
        super(SpiderConvert, self)._init( "__tmpl__" )
        self.output_files = [ self.map_file ]
    def _pre_exec( self ):
        self._make_script_file( 
            mrc_file=self.relpath( self.mrc_file )
        )



class SpiderDeleteFilledDensities( Spider ):
    args = [
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True }
    ]
    script_tmpl = "delete_filled_densities.spi"
    def _init( self, map_file, pdb_file, pixelsize, **kwargs ):
        self.map_file = self.abspath( map_file )
        self.pdb_file = self.abspath( pdb_file )
        self.pixelsize = pixelsize
        self.empty_map_file = self.outpath( "usermap.cpv" )
        super(SpiderDeleteFilledDensities, self)._init( "__tmpl__" )
        self.output_files = [ self.empty_map_file ]
    def _pre_exec( self ):
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            pdb_file=self.relpath( self.pdb_file ),
            pixelsize=self.pixelsize, 
            tmp_dir=self.relpath( self.output_dir ) + os.sep
        )



class SpiderBox( Spider ):
    args = [
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "res1", "type": "text", "help": "resno:chain, i.e. 10:A" },
        { "name": "res2", "type": "text" },
        { "name": "length", "type": "slider", "range": [1, 30] },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True },
        { "name": "resolution", "type": "slider", "range": [1, 10], "fixed": True, "help": "of the map_file" }
    ]
    script_tmpl = "box.spi"
    def _init( self, map_file, pdb_file, res1, res2, length, pixelsize, resolution, **kwargs ):
        self.map_file = self.abspath( map_file )
        self.pdb_file = self.abspath( pdb_file )
        self.res1 = res1
        self.res2 = res2
        self.length = length
        self.pixelsize = pixelsize
        self.resolution = resolution
        self.var_file = self.outpath( "variables.cpv" )
        self.box_file = self.outpath( "ergebnisse.cpv" )
        self.box_map_file = self.outpath( "boxil.cpv" )
        super(SpiderBox, self)._init( "__tmpl__" )
        self.output_files = [ self.box_file, self.box_map_file ]
    def _pre_exec( self ):
        coords1, coords2 = self._get_coords( self.pdb_file, self.res1, self.res2 )
        self._make_variables_file(
            coords1, coords2, self.length, self.pixelsize, self.resolution
        )
        self._make_script_file( 
            map_name=self.relpath( self.map_file, no_ext=True ), 
            var_name=self.relpath( self.var_file, no_ext=True )
        )
    def _get_coords( self, pdb_file, res1, res2 ):
        npdb = NumPdb( pdb_file )
        sele1 = numsele( res1 )
        sele2 = numsele( res2 )
        return npdb.center( **sele1 ), npdb.center( **sele2 )
    def _make_variables_file( self, coords1, coords2, length, pixelsize, resolution ):
        variables = "1 9 %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %i %4.2f %4.2f" % (
            coords1[0], coords1[1], coords1[2],
            coords2[0], coords2[1], coords2[2],
            length, pixelsize, resolution
        )
        with open( self.var_file, "w" ) as fp:
            fp.write( variables )


class SpiderReConvert( Spider ):
    args = [
        { "name": "box_file", "type": "file", "ext": "cpv" },
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "box_map_file", "type": "file", "ext": "cpv" },
    ]
    script_tmpl = "recon.spi"
    def _init( self, box_file, map_file, box_map_file, **kwargs ):
        self.box_file = self.abspath( box_file )
        self.map_file = self.abspath( map_file )
        self.box_map_file = self.abspath( box_map_file )
        self.mrc_file = self.outpath( "reconvert.mrc" )
        super(SpiderReConvert, self)._init( "__tmpl__" )
        self.output_files = [ self.mrc_file ]
    def _pre_exec( self ):
        self._make_script_file( 
            box_name=self.relpath( self.box_file, no_ext=True ),
            map_name=self.relpath( self.map_file, no_ext=True ),
            box_map_name=self.relpath( self.box_map_file, no_ext=True )
        )


class SpiderCrosscorrelation( Spider ):
    args = [
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "box_map_file", "type": "file", "ext": "cpv" },
        { "name": "box_file", "type": "file", "ext": "cpv" },
        { "name": "loop_file", "type": "file", "ext": "pdb" },
        { "name": "max_loops", "type": "slider", "range": [0, 200], "default_value": 100 }
    ]
    script_tmpl = "crosscorrelation.spi"
    def _init( self, map_file, box_map_file, box_file, loop_file, max_loops=100, **kwargs ):
        self.map_file = self.abspath( map_file )
        self.box_map_file = self.abspath( box_map_file )
        self.box_file = self.abspath( box_file )
        self.loop_file = self.abspath( loop_file )
        self.max_loops = False if not max_loops else int(max_loops)
        self.loop_dir = self.subdir( "loops" )
        self.crosscorrel_file = self.outpath( "crosscorrelation.cpv" )
        self.crosscorrel_json = self.outpath( "crosscorrelation.json" )
        super(SpiderCrosscorrelation, self)._init( "__tmpl__" )
        self.output_files = [ self.crosscorrel_file ]
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
        { "name": "mrc_file", "type": "file", "ext": "mrc" },
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "loop_file", "type": "file", "ext": "pdb" },
        { "name": "res1", "type": "text" },
        { "name": "res2", "type": "text" },
        { "name": "length", "type": "slider", "range": [1, 30] },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True },
        { "name": "resolution", "type": "slider", "range": [1, 10], "fixed": True },
        { "name": "max_loops", "type": "slider", "range": [0, 200], "default_value": 100 }
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, mrc_file, pdb_file, loop_file, 
               res1, res2, length, pixelsize, resolution, max_loops=100, **kwargs ):
        self.mrc_file = self.abspath( mrc_file )
        self.pdb_file = self.abspath( pdb_file )
        self.res1 = res1
        self.res2 = res2
        self.cropped_pdb = self.outpath( "cropped.pdb" )
        self.loop_file = self.abspath( loop_file )
        self.spider_convert = SpiderConvert( 
            self.mrc_file, **copy_dict( kwargs, run=False, output_dir=self.subdir("convert") )
        )
        self.spider_delete_filled_densities = SpiderDeleteFilledDensities( 
            self.spider_convert.map_file, self.cropped_pdb, pixelsize,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("delete_filled_densities") )
        )
        self.spider_box = SpiderBox(
            self.spider_delete_filled_densities.empty_map_file, 
            self.cropped_pdb, res1, res2, length, pixelsize, resolution,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("box") )
        )
        self.spider_reconvert = SpiderReConvert(
            self.spider_box.box_file,
            self.spider_convert.map_file,
            self.spider_box.box_map_file,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("reconvert") )
        )
        self.spider_crosscorrelation = SpiderCrosscorrelation(
            self.spider_convert.map_file, 
            self.spider_box.box_map_file, 
            self.spider_box.box_file, 
            self.loop_file,
            max_loops=max_loops,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("crosscorrelation") )
        )
        self.output_files = list( itertools.chain(
            self.spider_convert.output_files, 
            self.spider_delete_filled_densities.output_files,
            self.spider_box.output_files,
            self.spider_reconvert.output_files,
            self.spider_crosscorrelation.output_files,
            [ self.cropped_pdb ]
        ))
    def func( self ):
        self._crop_pdb()
        self.spider_convert()
        self.spider_delete_filled_densities()
        self.spider_box()
        self.spider_reconvert()
        self.spider_crosscorrelation()
    def _crop_pdb( self ):
        npdb = NumPdb( self.pdb_file )
        sele1 = numsele( self.res1 )
        sele2 = numsele( self.res2 )
        npdb.write( 
            self.cropped_pdb, 
            chain=sele1["chain"], 
            resno=[ sele1["resno"]+1, sele2["resno"]-1 ],
            invert=True
        )



