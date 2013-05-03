#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division



import re
import os
import shutil
import argparse

from utils.tool import CmdTool


PDB2PQR_CMD = "pdb2pqr.py"
APBS_CMD = "apbs"


# pdb2pqr.py --ff PARSE --apbs-input 1u19.pdb 1u19.pqr
# => 1u19-input.p 1u19.pqr 1u19.in

class Pdb2pqr( CmdTool ):
    def __init__( self, pdb_file, **kw ):
    	self.pdb_file = os.path.abspath( pdb_file )
    	stem = os.path.splitext( os.path.split( self.pdb_file )[-1] )[0]
    	print stem
        self.pqr_file = "%s.pqr" % stem
        self.apbsin_file = "%s.in" % stem
        self.apbsin_pickle = "%s-input.p" % stem
        self.cmd = [ 
            PDB2PQR_CMD, "--ff", "PARSE", "--apbs-input", self.pdb_file, self.pqr_file
        ]
        self.output_files = [ 
        	self.pqr_file, self.apbsin_file, self.apbsin_pickle
    	]
        super(Pdb2pqr, self).__init__( **kw )



# apbs 1u19.in
# => pot-PE0.dx io.mc

class Apbs( CmdTool ):
    def __init__( self, pdb_file, **kw ):
        self.pdb2pqr = Pdb2pqr( 
            pdb_file, output_dir=kw.get("output_dir"), 
            timeout=kw.get("timeout"), run=False
        )
        self.cmd = [ 
            APBS_CMD, self.pdb2pqr.apbsin_file
        ]
        self.output_files = self.pdb2pqr.output_files + \
            [ "pot-PE0", "io.mc" ]
        super(Apbs, self).__init__( **kw )
    def _pre_exec( self ):
        self.pdb2pqr()




def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-pdb', help='pdb file', type=str)
    parser.add_argument(
        '-output_dir', help='output dir', type=str, default=".")


    # parse the command line
    args = parser.parse_args()


    if args.pdb and args.output_dir:
    	Apbs( args.pdb, output_dir=args.output_dir )





if __name__ == "__main__":
    main()