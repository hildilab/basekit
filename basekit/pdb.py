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





def pdb_split( inputfile, output_dir ):
    """ author: Johanna Thiemann
        This function puts pdb-models into their own file.
    """ 
    BASE = re.sub("\..*$", "", re.sub("/.*$", "", inputfile[::-1])[::-1])
    with open( inputfile, "r" ) as fp_in:
        number=0
        for zeile in fp_in:
            if zeile.startswith('MODEL'):
                number=number+1
                newfile=open(outputDir+'/'+str(number).zfill(2)+'.pdb', 'w')
                newfilebb=open(outputDir+'/'+str(number).zfill(2)+'_bb.pdb', 'w')
            elif zeile.startswith('ATOM'):
                newfile.write(zeile)
                type=(' N ',' C ', ' CA ',' O ')
                for elem in type:
                    if elem in zeile: newfilebb.write(zeile)
            elif zeile.startswith('ENDMDL'):
                newfile.write('END')
                newfile.close()
                newfilebb.write('END')
                newfilebb.close()


class PdbSplit( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    def _init( self, pdb_file, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        self.output_files = [] # TODO doesn't work here
    def func( self ):
        pdb_split( self.pdb_file, self.output_dir )





def numpdb_test( pdb_file ):
    with Timer("read/parse pdb plain"):
        NumPdb( pdb_file, features={"phi_psi": False, "sstruc": False} )
    with Timer("read/parse pdb"):
        npdb = NumPdb( pdb_file )
    with Timer("dist"):
        print npdb.dist( {"chain":"A"}, {"chain":"B"} )
    with Timer("access phi/psi"):
        print np.nansum( npdb['phi'] )
    with Timer("sequence"):
        print npdb.sequence( chain="A" )
    with Timer("sstruc iter"):
        for numa in npdb.iter_sstruc():
            pass
    with Timer("resno iter"):
        for numa in npdb.iter_resno():
            pass


class NumpdbTest( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    no_output = True
    def _init( self, pdb_file, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        self.output_files = [] # TODO doesn't work here
    def func( self ):
        numpdb_test( self.pdb_file )

