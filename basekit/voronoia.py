from __future__ import with_statement
from __future__ import division



import os
import json

import utils.path
from utils import copy_dict
from utils.tool import CmdTool, ProviMixin

import provi_prep as provi



DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "voronoia" )
VOLUME_CMD = os.path.join( TMPL_DIR, "get_volume.exe" )


# get_volume.exe ex:0.1 rad:protor i:file.pdb o:out.vol

class Voronoia( CmdTool, ProviMixin ):
    """A wrapper around the 'voronoia' aka 'get_volume' programm."""
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "ex", "type": "slider", "range": [1, 5], "factor": 10, 
            "fixed": True, "default": 0.1  },
        { "name": "radii", "type": "select", "options": ["protor"], 
            "default": "protor" }
    ]
    out = [
        { "name": "vol_file", "file": "{pdb_file.stem}.vol" }
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "voronoia.provi"
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            "wine", VOLUME_CMD, 
            "ex:%0.1f"%float(self.ex), 
            "rad:%s"%self.radii,
            "i:%s"%self.pdb_file, 
            "o:%s"%self.vol_file
        ]
        self.output_files = [ 
        	self.vol_file
    	]
    def _post_exec( self ):
        provi.prep_volume( self.vol_file, self.pdb_file )
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            vol_file=self.relpath( self.vol_file )
        )



