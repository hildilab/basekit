from __future__ import with_statement
from __future__ import division



import os

import utils.path
from utils import copy_dict
from utils.tool import _, _dir_init, CmdTool

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "apbs" )

PDB2PQR_CMD = "pdb2pqr.py"
APBS_CMD = "apbs"


# pdb2pqr.py --ff PARSE --apbs-input 1u19.pdb 1u19.pqr
# => 1u19-input.p 1u19.pqr 1u19.in

class Pdb2pqr( CmdTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    out = [
        _( "pqr_file", file="{pdb_file.stem}.pqr" ),
        _( "apbsin_file", file="{pdb_file.stem}.in" ),
        _( "apbsin_pickle", file="{pdb_file.stem}-input.p" )
    ]
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            PDB2PQR_CMD, "--ff", "PARSE", "--apbs-input", 
            self.pdb_file, self.pqr_file
        ]


# apbs 1u19.in
# => pot-PE0.dx io.mc

class Apbs( CmdTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    out = [
        _( "dx_file", file="pot-PE0.dx" ),
        _( "mc_file", file="io.mc" )
    ]
    def _init( self, *args, **kwargs ):
        self.pdb2pqr = Pdb2pqr( 
            self.pdb_file, **copy_dict( kwargs, run=False ) 
        )
        self.cmd = [ 
            APBS_CMD,
            self.relpath( self.pdb2pqr.apbsin_file )
        ]
        self.output_files.extend( self.pdb2pqr.output_files )
    def _pre_exec( self ):
        self.pdb2pqr()



