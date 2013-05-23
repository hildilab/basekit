from __future__ import with_statement
from __future__ import division



import os

from utils import copy_dict
from utils.tool import CmdTool

import provi_prep as provi



DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
VOLUME_CMD = os.path.join( PARENT_DIR, "data", "voronoia", "get_volume.exe" )


# get_volume.exe ex:0.1 rad:protor i:file.pdb o:out.vol

class Voronoia( CmdTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "ex", "type": "slider", "range": [1, 5], "factor": 10, "fixed": True, "default_value": 0.1  },
        { "name": "radii", "type": "select", "options": ["protor"], "default_value": "protor" }

    ]
    def _init( self, pdb_file, ex=0.1, radii="protor", **kwargs ):
    	self.pdb_file = os.path.abspath( pdb_file )
    	stem = os.path.splitext( os.path.split( self.pdb_file )[-1] )[0]
        self.vol_file = "%s.vol" % stem
        self.cmd = [ 
            "wine", VOLUME_CMD, "ex:%0.1f"%float(ex), "rad:%s"%radii,
            "i:%s"%self.pdb_file, "o:%s"%self.vol_file
        ]
        self.output_files = [ 
        	self.vol_file
    	]
    def _post_exec( self ):
        provi.prep_volume(  self.vol_file, self.pdb_file )
