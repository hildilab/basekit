from __future__ import with_statement
from __future__ import division



import os
import re
import logging
import collections

from utils import memoize_m
from utils.tool import _, _dir_init, CmdTool, ProviMixin

import provi_prep as provi


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "voronoia" )
VOLUME_CMD = os.path.join( TMPL_DIR, "get_volume.exe" )

logging.basicConfig()
LOG = logging.getLogger('voronoia')


HOLE_PARTLY_FILLED = 1
HOLE_NOT_FILLED = 2
HOLE_PARTLY_FILLED_HETS_REMOVED = 3
HOLE_FILLED_HETS_REMOVED = 4

HoleNeighbour = collections.namedtuple( "HoleNeighbour", [
    "atomno", "atomname", "resno", "resname", "chain"
])
VolHole = collections.namedtuple( "VolHole", [
    "no", "type", "neighbours"
])

def parse_vol( vol_file, pdb_file ):
    pdb_coord_dict = provi.get_pdb_coord_dict( pdb_file )
    pdb_index_dict = provi.get_pdb_index_dict( pdb_file )
    vol_index_dict = {}

    holes = collections.OrderedDict()
    hole_types = []
    hole_list = []
    nrholes = {}
    
    i = 1
    with open( vol_file, "r" ) as fp:
        for l in fp:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                key = ( float(l[30:38]), float(l[38:46]), float(l[46:54]) )
                if key in pdb_coord_dict:
                    j = pdb_coord_dict[key]
                    vol_index_dict[ i ] = j
                    i += 1
                else:
                    raise Exception( 
                        "vol coords not in pdb coords dict. %s" % str(key) 
                    )
            elif l.startswith('HOLE NUMBER'):
                ls = map( int, l[12:].split() )
                holes[ ls[0] ] = ls[1:]
            elif l.startswith('NRHOLE'):
                p = ( ".* (\d+) \(holes partly filled\)"
                        ".* (\d+) \(holes not filled\)"
                        ".* (\d+) \(partly filled holes with Hets removed\)"
                        ".* (\d+) \(filled holes with Hets removed\).*" )
                m = re.match( p, l )
                if m:
                    nrholes = {
                        "partly_filled": int( m.group(1) ),
                        "not_filled": int( m.group(2) ),
                        "partly_filled_hets_removed": int( m.group(3) ),
                        "filled_hets_removed": int( m.group(4) ),
                    }
                    hole_types = ( 
                        [HOLE_PARTLY_FILLED] * int( m.group(1) ) +
                        [HOLE_NOT_FILLED] * int( m.group(2) ) +
                        [HOLE_PARTLY_FILLED_HETS_REMOVED] * int( m.group(3) ) +
                        [HOLE_FILLED_HETS_REMOVED] * int( m.group(4) )
                    )
                else:
                    raise Exception( "error parsing nrholes record" )

    if len(pdb_coord_dict) != len(vol_index_dict):
        LOG.debug( "number of atoms in vol and pdb coord dicts differs." )

    for no, h in holes.iteritems():
        neighbours = []
        for nb in h:
            if nb in vol_index_dict:
                neighbours.append( vol_index_dict[nb] )
            else:
                raise Exception( "hole neighbour index not found. %s" % nb )
        neighbours.sort()
        neighbours2 = []
        for nb in neighbours:
            pl = pdb_index_dict[ nb ]
            neighbours2.append(
                HoleNeighbour(
                    int( pl[6:11] ),    # atomno
                    pl[12:16],          # atomname
                    int( pl[22:26] ),   # resno
                    pl[17:20],          # resname
                    pl[21:22],          # chain
                )
            )
        hole_list.append( 
            VolHole( no, hole_types[ no-1 ], neighbours2 )
        )
        
    return {
        "nrholes": nrholes,
        "holes": hole_list
    }



# get_volume.exe ex:0.1 rad:protor i:file.pdb o:out.vol
class Voronoia( CmdTool, ProviMixin ):
    """A wrapper around the 'voronoia' aka 'get_volume' programm."""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "ex", type="float", range=[0.01, 0.5], step=0.01, default=0.1 ),
        _( "radii", type="select", options=["protor"], default="protor" )
    ]
    out = [
        _( "vol_file", file="{pdb_file.stem}.vol" ),
        _( "log_file", file="{pdb_file.stem}.log" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "voronoia.provi"
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            "wine", VOLUME_CMD, 
            "ex:%0.1f"%float(self.ex), 
            "rad:%s"%self.radii,
            "x:yes",
            "l:%s" %self.log_file,
            "i:%s"%self.pdb_file, 
            "o:%s"%self.vol_file
        ]
    def _post_exec( self ):
        provi.prep_volume( self.vol_file, self.pdb_file )
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            vol_file=self.relpath( self.vol_file )
        )
    @memoize_m
    def get_vol( self ):
        return parse_vol( self.vol_file, self.pdb_file )



