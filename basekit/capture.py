from __future__ import with_statement

import os
import urllib2
import collections

from poster.encode import multipart_encode
from poster.streaminghttp import register_openers

import numpy as np
np.seterr( all="raise" )

import utils.path
from utils.tool import (
    _, _dir_init, PyTool, RecordsMixin, ParallelMixin, ProviMixin
)



DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "capture" )
CAPTURE_URL = "http://capture.caltech.edu/"




def capture_web(pdb_file, output_file):
    """queries the CaPTURE webservice"""
    # Register the streaming http handlers with urllib2
    register_openers()
    # use relpath to hide local path
    with open( os.path.relpath( pdb_file ), "r" ) as fp:
        # headers contains the necessary Content-Type and Content-Length
        # datagen is a generator object that yields the encoded parameters
        datagen, headers = multipart_encode({
            "upfile": fp,
            "GO": "GO",
            "note": "note"
        })
        # Create the Request object
        request = urllib2.Request(CAPTURE_URL + "capture_ul.cgi", datagen, headers)
        # Actually do the request, get and read the response
        response = urllib2.urlopen(request).read()
    with open( output_file, 'w' ) as fp:
        fp.write( response )
    

def parse_capture( capture_file, pdb_id=None ):
    capture_list = []
    with open( capture_file, "r" ) as fp:
        for line in fp:
            if line[0:3] in ( "ARG", "LYS" ):
                x = line.strip().split("\t")
                capture_list.append([
                    x[0], int(x[1]), x[2],
                    x[3], int(x[4]), x[5],
                    float(x[6]), float(x[8])
                ])
    return capture_list


CaptureRecord = collections.namedtuple( 'CaptureRecord', [
    'pdb_id', 
    'resname_cation', 'resno_cation', 'chain_cation', 
    'resname_pi', 'resno_pi', 'chain_pi', 
    'es', "vdw"
])


class Capture( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """ CaPTURE web service wrapper (http://capture.caltech.edu/).
    """
    args = [
        _( "pdb_input", type="file", ext="pdb" ),
        _( "parse_only", type="checkbox", default=False )
    ]
    out = [
        _( "capture_file", file="CaPTURE_{pdb_input.stem}.txt" ),
    ]
    RecordsClass = CaptureRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "capture.provi"
    def _init( self, *args, **kwargs ):
        self.pdb_id = utils.path.stem( self.pdb_input )
        self._init_records( None, **kwargs )
        self._init_parallel( self.pdb_input, **kwargs )
    def func( self ):
        if not self.parse_only:
            capture_web( self.pdb_input, self.capture_file )
        self.records = [
            CaptureRecord( self.pdb_id, *x ) for x in
            parse_capture( self.capture_file )
        ]
        self.write()
    def _post_exec( self ):
        script = []
        for r in self.records:
            p = (
                 r.es, r.chain_cation, r.resno_cation,
                r.chain_pi, r.resno_pi
            )
            s = (
                "var colr = color( 'rwb', -5, 5, %f ); "
                "select { "
                    "(chain='%s' and resno=%i) or "
                    "(chain='%s' and resno=%i) "
                "}; "
                "wireframe 0.2; color @colr; " % p
            )
            script.append( s )
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_input ),
            script=" ".join( script )
        )

