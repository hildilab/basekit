from __future__ import with_statement
from __future__ import division

import os

from utils.tool import PyTool


def find_all( pdb_file, motif_type, output_dir ):
    print pdb_file, motif_type, output_dir


class CapsMotifFinder( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "motif_type", "type": "text", "default_value": "alpha_beta" }
    ]
    def _init( self, pdb_file, motif_type="alpha_beta", **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        self.motif_type = motif_type.split(",")
        self.output_files = [ "foo.txt" ]
    def func( self ):
        # self.output_dir is created by the base class
        find_all( self.pdb_file, self.motif_type, self.output_dir )
