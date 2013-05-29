from __future__ import with_statement
from __future__ import division

import os
import itertools
from string import Template

from utils import copy_dict
from utils.tool import Tool, CmdTool, ScriptMixin
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
    spi_tmpl_file = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, script_ext="spi", data_ext="cpv", **kwargs ):
        if script_file=="__tmpl__":
            script_file = os.path.join( self.output_dir, self.tmpl_file )
        self.script_file = os.path.abspath( script_file )
        print self.script_file
        # spider spi/cpv @box
        self.cmd = [
            SPIDER_CMD, 
            "%s/%s" % ( script_ext, data_ext ), 
            "@%s" % os.path.splitext( os.path.relpath( self.script_file, self.output_dir ) )[0]
        ]



class SpiderConvert( Spider ):
    args = [
        { "name": "mrc_file", "type": "file", "ext": "mrc" }
    ]
    tmpl_file = "convert.spi"
    def _init( self, mrc_file, **kwargs ):
        self.mrc_file = os.path.abspath( mrc_file )
        self.map_file = os.path.join( self.output_dir, "mapupload.cpv" )
        super(SpiderConvert, self)._init( "__tmpl__" )
        self.output_files = [ self.map_file ]
    def _pre_exec( self ):
        self._make_script_file( 
            mrc_file=os.path.relpath( self.mrc_file )
        )



class SpiderDeleteFilledDensities( Spider ):
    args = [
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True }
    ]
    tmpl_file = "delete_filled_densities.spi"
    def _init( self, map_file, pdb_file, pixelsize, **kwargs ):
        self.map_file = os.path.abspath( map_file )
        self.pdb_file = os.path.abspath( pdb_file )
        self.pixelsize = pixelsize
        self.empty_map_file = os.path.join( self.output_dir, "usermap.cpv" )
        super(SpiderDeleteFilledDensities, self)._init( "__tmpl__" )
        self.output_files = [ self.empty_map_file ]
    def _pre_exec( self ):
        self._make_script_file( 
            map_name=os.path.splitext( os.path.relpath( self.map_file ) )[0], 
            pdb_file=os.path.relpath( self.pdb_file ),
            pixelsize=self.pixelsize, 
            tmp_dir=os.path.relpath( self.output_dir ) + os.sep
        )


class SpiderBox( Spider ):
    args = [
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "res1", "type": "text" },
        { "name": "res2", "type": "text" },
        { "name": "length", "type": "slider", "range": [1, 30] },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True },
        { "name": "resolution", "type": "slider", "range": [1, 10], "fixed": True }
    ]
    tmpl_file = "box.spi"
    def _init( self, map_file, pdb_file, res1, res2, length, pixelsize, resolution, **kwargs ):
        self.map_file = os.path.abspath( map_file )
        self.pdb_file = os.path.abspath( pdb_file )
        self.res1 = res1
        self.res2 = res2
        self.length = length
        self.pixelsize = pixelsize
        self.resolution = resolution
        self.var_file = os.path.join( self.output_dir, "variables.cpv" )
        self.box_file = os.path.join( self.output_dir, "ergebnisse.cpv" )
        self.box_map_file = os.path.join( self.output_dir, "boxil.cpv" )
        super(SpiderBox, self)._init( "__tmpl__" )
        self.output_files = [ self.box_file, self.box_map_file ]
    def _pre_exec( self ):
        coords1, coords2 = self._get_coords( self.pdb_file, self.res1, self.res2 )
        self._make_variables_file(
            coords1, coords2, self.length, self.pixelsize, self.resolution
        )
        self._make_script_file( 
            map_name=os.path.splitext( os.path.relpath( self.map_file ) )[0], 
            var_name=os.path.splitext( os.path.relpath( self.var_file ) )[0]
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



class SpiderCrosscorrelation( Spider ):
    args = [
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "box_map_file", "type": "file", "ext": "cpv" },
        { "name": "box_file", "type": "file", "ext": "cpv" },
        { "name": "loop_file", "type": "file", "ext": "pdb" },
        { "name": "max_loops", "type": "slider", "range": [0, 200], "default_value": 100 }
    ]
    tmpl_file = "crosscorrelation.spi"
    def _init( self, map_file, box_map_file, box_file, loop_file, max_loops=100, **kwargs ):
        self.map_file = os.path.abspath( map_file )
        self.box_map_file = os.path.abspath( box_map_file )
        self.box_file = os.path.abspath( box_file )
        self.loop_file = os.path.abspath( loop_file )
        self.max_loops = False if not max_loops else int(max_loops)
        self.loop_dir = os.path.join( self.output_dir, "loops" )
        self.crosscorrel_file = os.path.join( self.output_dir, "crosscorrelation.cpv" )
        super(SpiderCrosscorrelation, self)._init( "__tmpl__" )
        self.output_files = [ self.crosscorrel_file ]
    def _pre_exec( self ):
        self._split_loop_file()
        self._make_script_file( 
            map_name=os.path.splitext( os.path.relpath( self.map_file ) )[0], 
            box_map_name=os.path.splitext( os.path.relpath( self.box_map_file ) )[0],
            box_name=os.path.splitext( os.path.relpath( self.box_file ) )[0],
            loop_dir=os.path.relpath( self.loop_dir ) + os.sep,
            max_loops=self.max_loops or 999
        )
    def _split_loop_file( self ):
        PdbSplit( 
            self.loop_file, output_dir=self.loop_dir, backbone_only=True, 
            max_models=self.max_loops, resno_ignore=[ 1000, 2000 ], zfill=3
        )


class LoopCrosscorrel( Tool ):
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
    def _init( self, mrc_file, pdb_file, loop_file, 
               res1, res2, length, pixelsize, resolution, max_loops=100, **kwargs ):
        def outdir( subdir ):
            return os.path.join( self.output_dir, subdir )
        self.mrc_file = os.path.abspath( mrc_file )
        self.pdb_file = os.path.abspath( pdb_file )
        self.loop_file = os.path.abspath( loop_file )
        self.spider_convert = SpiderConvert( 
            self.mrc_file, **copy_dict( kwargs, run=False, output_dir=outdir("convert") )
        )
        self.spider_delete_filled_densities = SpiderDeleteFilledDensities( 
            self.spider_convert.map_file, self.pdb_file, pixelsize,
            **copy_dict( kwargs, run=False, output_dir=outdir("delete_filled_densities") )
        )
        self.spider_box = SpiderBox(
            self.spider_delete_filled_densities.empty_map_file, 
            self.pdb_file, res1, res2, length, pixelsize, resolution,
            **copy_dict( kwargs, run=False, output_dir=outdir("box") )
        )
        self.spider_crosscorrelation = SpiderCrosscorrelation(
            self.spider_convert.map_file, 
            self.spider_box.box_map_file, 
            self.spider_box.box_file, self.loop_file,
            max_loops=max_loops,
            **copy_dict( kwargs, run=False, output_dir=outdir("crosscorrelation") )   
        )
        self.output_files = list( itertools.chain(
            self.spider_convert.output_files, 
            self.spider_delete_filled_densities.output_files,
            self.spider_box.output_files, 
            self.spider_crosscorrelation.output_files
        ))
    def _pre_exec( self ):
        self.spider_convert()
        self.spider_delete_filled_densities()
        self.spider_box()
        self.spider_crosscorrelation()



