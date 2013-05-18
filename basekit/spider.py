from __future__ import with_statement
from __future__ import division

import os
from string import Template

from utils.tool import CmdTool, ScriptMixin, make_args
from utils.numpdb import NumPdb, numsele


DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "spider" )
SPIDER_CMD = "spider" 


def fabs( fpath ):
    return os.path.abspath( fpath )

def fname( fpath ):
    return os.path.splitext( fabs( fpath ) )[0]



class Spider( CmdTool, ScriptMixin ):
    args = make_args([
        { "name": "script_file", "type": "file", "ext": "spi" }
    ])
    spi_tmpl_file = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, script_ext="spi", data_ext="cpv", **kwargs ):
        # spider spi/cpv @box
        self.script_file = fabs( script_file )
        self.cmd = [
            SPIDER_CMD, "%s/%s" % ( script_ext, data_ext ), "@%s" % fname( script_file )
        ]



class SpiderConvert( Spider ):
    args = make_args([
        { "name": "mrc_file", "type": "file", "ext": "mrc" }
    ])
    tmpl_file = "convert.spi"
    def _init( self, mrc_file, **kwargs ):
        script_file = self._make_script_file( mrc_file=fabs( mrc_file ) )
        super(SpiderConvert, self)._init( script_file )
        self.output_files = [ "mapupload.cpv" ]



class SpiderDeleteFilledDensities( Spider ):
    args = make_args([
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True }
    ])
    tmpl_file = "delete_filled_densities.spi"
    def _init( self, map_file, pdb_file, pixelsize=1, **kwargs ):
        script_file = self._make_script_file( 
            map_name=fname( map_file ), pdb_file=fabs( pdb_file ), pixelsize=pixelsize,
            tmp_dir=self.output_dir
        )
        super(SpiderDeleteFilledDensities, self)._init( script_file )
        self.output_files = [ "usermap.cpv" ]



class SpiderBox( Spider ):
    args = make_args([
        { "name": "map_file", "type": "file", "ext": "cpv" },
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "res1", "type": "text" },
        { "name": "res2", "type": "text" },
        { "name": "length", "type": "slider", "range": [1, 30] },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True },
        { "name": "resolution", "type": "slider", "range": [1, 10], "fixed": True }
    ])
    tmpl_file = "box.spi"
    def _init( self, map_file, pdb_file, res1, res2, length, pixelsize, resolution, **kwargs ):
        coords1, coords2 = self._get_coords( fabs( pdb_file ), res1, res2 )
        var_file = self._make_variables_file(
            coords1, coords2, length, pixelsize, resolution
        )
        script_file = self._make_script_file( 
            map_name=fname( map_file ), var_name=fname( var_file )
        )
        super(SpiderBox, self)._init( script_file )
        self.output_files = [ "ergebnisse.cpv" ]
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
        var_file = os.path.join( self.output_dir, "variables.cpv" )
        with open( var_file, "w" ) as fp:
            fp.write( variables )
        return var_file



class SpiderCrosscorrelation( Spider ):
    args = make_args([
        { "name": "map_file1", "type": "file", "ext": "cpv" },
        { "name": "map_file2", "type": "file", "ext": "cpv" },
        { "name": "box_file", "type": "file", "ext": "cpv" },
        { "name": "loop_dir", "type": "text" }
    ])
    tmpl_file = "crosscorrelation.spi"
    def _init( self, map_file1, map_file2, box_file, loop_dir, **kwargs ):
        script_file = self._make_script_file( 
            map_name1=fname( map_file1 ), map_name2=fname( map_file2 ),
            box_name=fname( box_file ), loop_dir=fabs( loop_dir )+os.sep
        )
        super(SpiderCrosscorrelation, self)._init( script_file )
        self.output_files = [ "crosscorrelation.cpv" ]


