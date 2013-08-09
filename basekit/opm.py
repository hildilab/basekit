from __future__ import with_statement
from __future__ import division


import os
import json
import urllib2

import numpy as np


from utils import memoize_m
from utils.math import norm
from utils.tool import _, _dir_init, PyTool, ProviMixin



DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "opm" )
OPM_PDB_URL = "http://opm.phar.umich.edu/pdb/{pdb_id:}.pdb"
OPM_LOCAL_PATH = os.environ.get("OPM_LOCAL_PATH", "")


def opm( pdb_id ):
    pdb_id = pdb_id.lower()
    try:
        path = os.path.join( OPM_LOCAL_PATH, "%s.pdb" % pdb_id )
        with open( path, "r" ) as fp:
            return fp.read()
    except IOError:
        try:
            url = OPM_PDB_URL.format( pdb_id=pdb_id )
            return urllib2.urlopen( url ).read()
        except urllib2.HTTPError:
            raise Exception("Opm url error")


def ppm( pdb_file ):
    raise NotImplementedError


class OpmMixin( object ):
    def _post_exec( self ):
        with open( self.processed_file, "w" ) as fp:
            with open( self.opm_file, "r" ) as fp_opm:
                for line in fp_opm:
                    if line[17:20]!="DUM" and line[0:6]!="REMARK":
                        fp.write( line )
        self.make_mplane_file()
        self._make_provi_file(
            pdb_file=self.relpath( self.processed_file ),
            mplane_file=self.relpath( self.mplane_file )
        )
    @memoize_m
    def get_planes( self ):
        # TODO works only for two planes
        with open( self.opm_file ) as fp:
            coords={'O':[],'N':[]}
            for line in fp:
                if line[0:6]=="HETATM" and line[17:20]=="DUM":
                    atm = line[13:14]
                    c = np.array( map( float, [ 
                        line[30:38], line[38:46], line[46:54] 
                    ]))
                    if len(coords[ atm ])==2:
                        vn1 = norm( coords[ atm ][1] - coords[ atm ][0] )
                        vn2 = norm( c - coords[ atm ][0] )
                        if not np.allclose( vn1, vn2 ):
                            coords[ atm ].append(c)
                    elif len(coords[ atm ])<2:
                        coords[ atm ].append(c)
                if sum( map( len, coords.values() ) )==6:
                    break
            else:
                raise Exception( "could not find plane coordinates" )
        return np.array([ coords["N"], coords["O"] ])
    def make_mplane_file( self ):
        with open( self.mplane_file, "w" ) as fp:
            mp = self.get_planes().tolist()
            json.dump( mp, fp )


OPM_OUT = [
    _( "opm_file", file="{pdb_file.stem}_opm.pdb" ),
    _( "processed_file", file="{pdb_file.stem}_proc.pdb" ),
]


class Opm( OpmMixin, PyTool, ProviMixin ):
    """A tool to access the OPM database"""
    args = [
        _( "pdb_id", type="str" ),
    ]
    out = [
        _( "opm_file", file="{pdb_id}_opm.pdb" ),
        _( "mplane_file", file="{pdb_id}.mplane" ),
        _( "processed_file", file="{pdb_id}_proc.pdb" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "opm.provi"
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        with open( self.opm_file, "w" ) as fp:
            fp.write( opm( self.pdb_id ) )



class Ppm( OpmMixin, PyTool, ProviMixin ):
    """A tool to query the PPM webservice"""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
    ]
    out = [
        _( "opm_file", file="{pdb_file.stem}_opm.pdb" ),
        _( "mplane_file", file="{pdb_file.stem}.mplane" ),
        _( "processed_file", file="{pdb_file.stem}_proc.pdb" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "opm.provi"
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        ppm( self.pdb_file )

