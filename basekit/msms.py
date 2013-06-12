#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division


import re
import sys
import os

import utils.path
from utils import copy_dict
from utils.tool import CmdTool


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



class Msms( CmdTool ):
    """A wrapper around the MSMS program."""
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "density", "type": "slider", "range": [1, 10], 
            "fixed": True, "default": 1.0  }
    ]
    out = [
        { "name": "area_file", "file": "area.area" },
        { "name": "face_file", "file": "tri_surface.face" },
        { "name": "vert_file", "file": "tri_surface.vert" }
    ]
    def _init( self, *args, **kwargs ):
        self.pdb2xyzr = Pdb2xyzr( 
            self.pdb_file, **copy_dict( kwargs, run=False ) 
        )
        self.cmd = [ 
            MSMS_CMD, "-if", self.pdb2xyzr.xyzr_file, 
            "-af", "area", "-of", "tri_surface", 
            "-density", str(self.density)
        ]
        self.output_files.extend( self.pdb2xyzr.output_files )
    def _pre_exec( self ):
        self.pdb2xyzr()



class Pdb2xyzr( CmdTool ):
    """A pdb to xyzr format converter based on OpenBabel."""
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    out = [
        { "name": "pdb_prep_file", "file": "{pdb_file.stem}_prep.pdb" },
        { "name": "xyzr_file", "file": "{pdb_file.stem}.vol" }
    ]
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            BABEL_CMD, '-i', 'pdb', self.pdb_prep_file,
            '-o', 'msms', self.xyzr_file 
        ]
        self.output_files = [ self.pdb_prep_file, self.xyzr_file ]
    def _pre_exec( self ):
        pdb_select( self.pdb_file, self.pdb_prep_file )


