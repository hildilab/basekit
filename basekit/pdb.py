from __future__ import with_statement

import os
import re
import urllib2

from utils.tool import PyTool, make_args




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
    args = make_args([
        { "name": "pdb_id", "type": "text" }
    ])
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
    args = make_args([
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ])
    def _init( self, pdb_file, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        self.output_files = [] # TODO doesn't work here
    def func( self ):
        pdb_split( self.pdb_file, self.output_dir )

