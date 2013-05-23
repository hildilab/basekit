from __future__ import with_statement
from __future__ import division



import os

from utils import copy_dict
from utils.tool import CmdTool


PDB2PQR_CMD = "pdb2pqr.py"
APBS_CMD = "apbs"


# pdb2pqr.py --ff PARSE --apbs-input 1u19.pdb 1u19.pqr
# => 1u19-input.p 1u19.pqr 1u19.in

class Pdb2pqr( CmdTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    def _init( self, pdb_file, **kwargs ):
    	self.pdb_file = os.path.abspath( pdb_file )
    	stem = os.path.splitext( os.path.split( self.pdb_file )[-1] )[0]
        self.pqr_file = "%s.pqr" % stem
        self.apbsin_file = "%s.in" % stem
        self.apbsin_pickle = "%s-input.p" % stem
        self.cmd = [ 
            PDB2PQR_CMD, "--ff", "PARSE", "--apbs-input", self.pdb_file, self.pqr_file
        ]
        self.output_files = [ 
        	self.pqr_file, self.apbsin_file, self.apbsin_pickle
    	]


# apbs 1u19.in
# => pot-PE0.dx io.mc

class Apbs( CmdTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    def _init( self, pdb_file, **kwargs ):
        self.pdb2pqr = Pdb2pqr( pdb_file, **copy_dict( kwargs, run=False ) )
        self.cmd = [ APBS_CMD, self.pdb2pqr.apbsin_file ]
        self.output_files = self.pdb2pqr.output_files + [ "pot-PE0.dx", "io.mc" ]
    def _pre_exec( self ):
        self.pdb2pqr()



