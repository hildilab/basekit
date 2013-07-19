#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division


import re
import sys
import os

import utils.path
from utils import copy_dict
from utils.tool import _, _dir_init, CmdTool, ProviMixin

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "msms" )

MSMS_CMD = "msms"
BABEL_CMD = "babel"



def pdb_select( input_pdb, output_pdb ):
    with open( input_pdb, "r" ) as fp_in:
        with open( output_pdb, "w" ) as fp_out:
            for line in fp_in:
                if line.startswith("ATOM"):
                    fp_out.write( line )
                elif line.startswith("ENDMDL"):
                    break



class Msms( CmdTool, ProviMixin ):
    """A wrapper around the MSMS program."""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "density", type="float", range=[0.5, 10], step=0.5, default=1.0 ),
        _( "hdensity", type="float", range=[1.0, 20], step=1.0, default=3.0 ),
        _( "all_components", type="checkbox", default=False ),
        _( "no_area", type="checkbox", default=False ),
    ]
    out = [
        _( "area_file", file="area.area" ),
        _( "face_file", file="tri_surface.face" ),
        _( "vert_file", file="tri_surface.vert" ),
        _( "provi_file", file="msms.provi" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "msms.provi"
    def _init( self, *args, **kwargs ):
        self.pdb2xyzr = Pdb2xyzr( 
            self.pdb_file, **copy_dict( kwargs, run=False ) 
        )
        self.cmd = [ 
            MSMS_CMD, "-if", self.pdb2xyzr.xyzr_file, 
            "-af", "area", "-of", "tri_surface", 
            "-density", str(self.density),
            "-hdensity", str(self.hdensity),
        ]
        if self.all_components:
            self.cmd.append( "-all_components" )
        if self.no_area:
            self.cmd.append( "-no_area" )
        self.output_files.extend( self.pdb2xyzr.output_files )
    def _pre_exec( self ):
        self.pdb2xyzr()
    def _post_exec( self ):
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            vert_file=self.relpath( self.vert_file )
        )



class Pdb2xyzr( CmdTool ):
    """A pdb to xyzr format converter based on OpenBabel."""
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    out = [
        _( "pdb_prep_file", file="{pdb_file.stem}_prep.pdb" ),
        _( "xyzr_file", file="{pdb_file.stem}.xyzr" )
    ]
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            BABEL_CMD, '-i', 'pdb', self.pdb_prep_file,
            '-o', 'msms', self.xyzr_file 
        ]
        self.output_files = [ self.pdb_prep_file, self.xyzr_file ]
    def _pre_exec( self ):
        pdb_select( self.pdb_file, self.pdb_prep_file )


