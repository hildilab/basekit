from __future__ import with_statement
from __future__ import division

import os

from utils.tool import PyTool, make_args


def find_all( tool ):
    print tool.output_dir, tool.pdb_file
    print tool.motif_type


class CapsMotifFinder( PyTool ):
    args = make_args([
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "motif_type", "type": "text", "default_value": "alpha_beta" }
    ])
    def _init( self, pdb_file, motif_type="alpha_beta", **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        self.motif_type = motif_type.split(",")
        self.func = find_all
        self.output_files = [ "foo.txt" ]
