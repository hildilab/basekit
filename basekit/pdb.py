from __future__ import with_statement

import os
import re
import urllib2
import gzip
import collections
import json

import numpy as np
np.seterr( all="raise" )

import utils.path
from utils.tool import PyTool
from utils.timer import Timer
from utils.numpdb import NumPdb, numdist, ATOMS
from utils.db import get_pdb_files



# RESIDUE   TPO     22
# CONECT      CG2    4 CB   HG21 HG22 HG23
def parse_het_dictionary( het_file ):
    het_dict = collections.defaultdict( list )
    key = None
    with open( het_file, "r" ) as fp:
        for line in fp:
            if line.startswith("RESIDUE"):
                key = line[10:13]
            if line.startswith("CONECT"):
                het_dict[ key ].append( line[11:15] )
    print len(het_dict)
    aminoacid_list = []
    for key, atom_list in het_dict.iteritems():
        if ATOMS['backbone'].issubset( atom_list ):
            aminoacid_list.append( key )
    print len(aminoacid_list)
    return aminoacid_list
    

class PdbHetDictionary( PyTool ):
    """
    Tool to extract the names of (non-standard) amino acids from
    ftp://ftp.rcsb.org/pub/pdb/data/monomers/het_dictionary.txt
    """
    args = [
        { "name": "het_file", "type": "file", "ext": "txt", "help": "hetero dictionary" }
    ]
    def _init( self, het_file, **kwargs ):
        self.het_file = self.abspath( het_file )
        self.aa_file = self.outpath( "aa.json" )
    def func( self ):
        aminoacid_list = parse_het_dictionary( self.het_file )
        with open( self.aa_file, "w" ) as fp:
            json.dump( aminoacid_list, fp )




def unzip_pdb( fpath ):
    fdir, fname = os.path.split( fpath )
    fpdb = os.path.join( fdir, fname[3:7]+".pdb" )
    with open( fpdb, "w" ) as fp:
        with gzip.open( fpath, "rb" ) as zp:
            fp.write( zp.read() )
    return True

class PdbUnzip( PyTool ):
    args = [
        { "name": "pdb_archive", "type": "text" }
    ]
    no_output = True
    def _init( self, pdb_archive, **kwargs ):
        self.pdb_archive = self.abspath( pdb_archive )
    def func( self ):
        pdb_file_list = get_pdb_files( self.pdb_archive, pattern="ent.gz" )
        for pdb_file in pdb_file_list:
            unzip_pdb( pdb_file )




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
        { "name": "pdb_id", "type": "text", 
          "help": "single id or multiple, seperated by spaces or commas" }
    ]
    def _init( self, pdb_id, **kwargs ):
        self.pdb_id_list = map( lambda x: x[0:4], re.split("[\s,]+", pdb_id.strip() ) )
        self.pdb_file_list = map(
            lambda x: self.outpath( "%s.pdb" % x ),
            self.pdb_id_list
        )    
        self.output_files = [] + self.pdb_file_list
    def func( self ):
        for pdb_id, pdb_file in zip( self.pdb_id_list, self.pdb_file_list ):
            pdb_download( pdb_id, pdb_file )





def pdb_split( pdb_file, output_dir, backbone_only=False, 
               resno_ignore=False, max_models=False, zfill=False ):
    """ author: Johanna Thiemann
        author: Alexander Rose
        This function puts pdb-models into their own file.
    """ 
    backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
    bb_tag = "_bb" if backbone_only else ""
    model_no = 0
    if max_models and not zfill:
        zfill = len( str(max_models) )
    with open( pdb_file, "r" ) as fp:
        for line in fp:
            if line[0:5]=='MODEL':
                model_no += 1
                if max_models and model_no>max_models:
                    break
                file_name = "%s%s.pdb" % ( str(model_no).zfill( zfill ), bb_tag )
                file_path = os.path.join( output_dir, file_name )
                with open( file_path, 'w') as fp_out:
                    for line in fp:
                        if line[0:4]!='ATOM':
                            if line[0:6]=='ENDMDL':
                                fp_out.write( 'END' )
                                break
                            continue
                        if backbone_only and line[12:16] not in backbone:
                            continue
                        if resno_ignore and int( line[22:26] ) in resno_ignore:
                            continue
                        fp_out.write( line )
                        


class PdbSplit( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "backbone_only", "type": "checkbox", "default_value": False },
        { "name": "max_models", "type": "slider", "range": [0, 100], "default_value": 0 },
        { "name": "resno_ignore", "type": "text", "default_value": "" }
    ]
    def _init( self, pdb_file, backbone_only=False, max_models=False, 
               resno_ignore=False, zfill=False, **kwargs ):
        self.pdb_file = self.abspath( pdb_file )
        self.backbone_only = backbone_only
        self.max_models = int( max_models )
        self.zfill = int( zfill )
        self.resno_ignore = False
        if resno_ignore:
            if isinstance( resno_ignore, basestring ):
                self.resno_ignore = map( int, resno_ignore.split(",") )
            else:
                self.resno_ignore = resno_ignore
        self.output_files = [] # TODO doesn't work here
    def func( self ):
        pdb_split( 
            self.pdb_file, self.output_dir, 
            backbone_only=self.backbone_only,
            max_models=self.max_models,
            resno_ignore=self.resno_ignore,
            zfill=self.zfill
        )





def numpdb_test( pdb_file ):
    with Timer("read/parse pdb plain"):
        NumPdb( pdb_file, features={
            "phi_psi": False, 
            "sstruc": False, 
            "backbone_only": True,
            "detect_incomplete": False
        })
    with Timer("read/parse pdb incomplete"):
        NumPdb( pdb_file, features={"phi_psi": False, "sstruc": False, "backbone_only": False, "detect_incomplete": False} )
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
    with Timer("access phi/psi, altloc"):
        print np.nansum( npdb.get('phi', altloc=[" ", "A"] ) )
    with Timer("sstruc iter"):
        for numa in npdb.iter_sstruc():
            pass
    with Timer("resno iter"):
        for numa in npdb.iter_resno():
            pass
    with Timer("plain resno iter"):
        for numa in npdb._iter_resno():
            pass
    # with Timer("sequence"):
    #     first_chain = npdb["chain"][0]
    #     print npdb.sequence( chain=first_chain )
    with Timer("write pdb"):
        npdb.write( "test.pdb" )


class NumpdbTest( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    no_output = True
    def _init( self, pdb_file, **kwargs ):
        self.pdb_file = self.abspath( pdb_file )
    def func( self ):
        numpdb_test( self.pdb_file )

