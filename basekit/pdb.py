from __future__ import with_statement

import os
import re
import urllib2

import numpy as np

from utils.tool import PyTool
from utils.timer import Timer
from utils.numpdb import NumPdb, numdist




def pdb_download( pdb_id, output_file ):
    pdb_download_urls = [
        "http://www.rcsb.org/pdb/files/",
        "http://198.202.122.52/pdb/files/", # outdated?
        "http://198.202.122.51/pdb/files/"  # outdated?
    ]
    for base_url in pdb_download_urls:
        try:
            url = "%s%s.pdb" % ( base_url, pdb_id )
            request = urllib2.Request( url )
            response = urllib2.urlopen( request )
            with open( output_file, 'wb' ) as fp:
                fp.write( response.read() )
            return
        except:
            pass


class PdbDownload( PyTool ):
    args = [
        { "name": "pdb_id", "type": "text" }
    ]
    def _init( self, pdb_id, **kwargs ):
        self.pdb_id = pdb_id.strip()[0:4]
        self.pdb_file = os.path.join( self.output_dir, "%s.pdb" % self.pdb_id )
        self.output_files = [ self.pdb_file ]
    def func( self ):
        pdb_download( self.pdb_id, self.pdb_file )





def pdb_split( pdb_file, output_dir, backbone_only=False ):
    """ author: Johanna Thiemann
        author: Alexander Rose
        This function puts pdb-models into their own file.
    """ 
    backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
    bb_tag = "_bb" if backbone_only else ""
    model_no=0
    with open( pdb_file, "r" ) as fp:
        for line in fp:
            if line[0:5]=='MODEL':
                model_no += 1
                file_name = "%05i%s.pdb" % ( model_no, bb_tag )
                file_path = os.path.join( output_dir, file_name )
                with open( file_path, 'w') as fp_out:
                    for line in fp:
                        if line[0:4]=='ATOM' and ( not backbone_only or line[12:16] in backbone ):
                            fp_out.write( line )
                        elif line[0:6]=='ENDMDL':
                            fp_out.write( 'END' )
                            break


class PdbSplit( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "backbone_only", "type": "checkbox", "default_value": False }
    ]
    def _init( self, pdb_file, backbone_only=False, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        self.backbone_only = backbone_only
        self.output_files = [] # TODO doesn't work here
    def func( self ):
        pdb_split( self.pdb_file, self.output_dir, backbone_only=self.backbone_only )





def numpdb_test( pdb_file ):
    with Timer("read/parse pdb plain"):
        NumPdb( pdb_file, features={
            "phi_psi": False, 
            "sstruc": False, 
            "backbone_only": True,
            "detect_incomplete": False
        })
    with Timer("read/parse pdb incomplete"):
        NumPdb( pdb_file, features={"phi_psi": False, "sstruc": False, "backbone_only": True} )
    with Timer("read/parse pdb phi/psi"):
        NumPdb( pdb_file, features={"sstruc": False} )
    with Timer("read/parse pdb"):
        npdb = NumPdb( pdb_file )
    # with Timer("resno iter2"):
    #     for numa_list in npdb.iter_resno2( 6, chain="A", resno=[1,22] ):
    #         print [ numa["resno"][0] for numa in numa_list ]
    # with Timer("dist"):
    #     print npdb.dist( {"chain":"A"}, {"chain":"B"} )
    with Timer("access phi/psi"):
        print np.nansum( npdb['phi'] )
    # with Timer("sstruc iter"):
    #     for numa in npdb.iter_sstruc():
    #         pass
    # with Timer("resno iter"):
    #     for numa in npdb.iter_resno( chain="A" ):
    #         pass
    # with Timer("sequence"):
    #     first_chain = npdb["chain"][0]
    #     print npdb.sequence( chain=first_chain )


class NumpdbTest( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    no_output = True
    def _init( self, pdb_file, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
    def func( self ):
        numpdb_test( self.pdb_file )

