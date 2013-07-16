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
import utils.numpdb as numpdb
from utils import try_int
from utils.tool import _, _dir_init, PyTool, ProviMixin
from utils.timer import Timer
from utils.db import get_pdb_files


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "pdb" )



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
        if numpdb.ATOMS['backbone'].issubset( atom_list ):
            aminoacid_list.append( key )
    print len(aminoacid_list)
    return aminoacid_list
    
class PdbHetDictionary( PyTool ):
    """
    Tool to extract the names of (non-standard) amino acids from
    ftp://ftp.rcsb.org/pub/pdb/data/monomers/het_dictionary.txt
    """
    args = [
        _( "het_file", type="file", ext="txt", help="hetero dictionary" )
    ]
    out = [
        _( "aa_file", file="aa.json" )
    ]
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
        _( "pdb_archive", type="dir" )
    ]
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
        _( "pdb_id", type="text", 
            help="single id or multiple, seperated by spaces or commas" )
    ]
    def _init( self, *args, **kwargs ):
        self.pdb_id_list = map( 
            lambda x: x[0:4], 
            re.split("[\s,]+", self.pdb_id.strip() ) 
        )
        self.pdb_file_list = map(
            lambda x: self.outpath( "%s.pdb" % x ),
            self.pdb_id_list
        )    
        self.output_files.extend( self.pdb_file_list )
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
                        if resno_ignore:
                            if try_int( line[22:26] ) in resno_ignore:
                                continue
                        fp_out.write( line )
                        
class PdbSplit( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "backbone_only", type="checkbox", default=False ),
        _( "max_models", type="slider", range=[0, 100], default=0 ),
        _( "resno_ignore", type="text", default="" ),
        _( "zfill", type="slider", range=[0, 8], default=0 )
    ]
    def _init( self, *args, **kwargs ):
        if self.resno_ignore:
            if isinstance( self.resno_ignore, basestring ):
                self.resno_ignore = map( int, self.resno_ignore.split(",") )
    def func( self ):
        pdb_split( 
            self.pdb_file, self.output_dir, 
            backbone_only=self.backbone_only,
            max_models=self.max_models,
            resno_ignore=self.resno_ignore,
            zfill=self.zfill
        )



class PdbEdit( PyTool ):
    """ Edits a pdb file. Manipulations are done in this order:
        center, shift, box.
    """
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "center", type="checkbox", default=False ),
        _( "shift", type="float", nargs=3, default=None,
            metavar=("X", "Y", "Z") ),
        _( "box", type="float", nargs=6, default=None,
            help="Two sets of 3d coordinates; first one corner, "
                "then the extent at which the opposing corner is. "
                "Outputs only the atoms within the box.",
            metavar=("X", "Y", "Z", "EX", "EY", "EZ") ),
    ]
    out = [
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    def func( self ):
        npdb = numpdb.NumPdb( self.pdb_file, {
            "phi_psi": False,
            "sstruc": False,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False
        })
        sele = None

        if self.center:
            npdb['xyz'] -= npdb['xyz'].mean( axis=0 )

        if self.shift:
            shift = np.array( self.shift[0:3] )
            npdb['xyz'] += shift

        if self.box:
            corner1 = np.array( self.box[0:3] )
            corner2 = corner1 + np.array( self.box[3:6] )
            sele = (
                ( npdb['xyz']>corner1 ) & 
                ( npdb['xyz']<corner2 )
            ).all( axis=1 )

        npdb.copy( sele=sele ).write( self.edited_pdb_file )




class PdbSuperpose( PyTool, ProviMixin ):
    args = [
        _( "pdb_file1", type="file", ext="pdb" ),
        _( "pdb_file2", type="file", ext="pdb" ),
        _( "sele1", type="sele" ),
        _( "sele2", type="sele" ),
        _( "subset", type="text", default="CA" )
    ]
    out = [
        _( "superposed_file", file="superposed.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "superpose.provi"
    def func( self ):
        npdb1 = numpdb.NumPdb( self.pdb_file1 )
        npdb2 = numpdb.NumPdb( self.pdb_file2 )
        numpdb.superpose( 
            npdb1, npdb2, self.sele1, self.sele2, 
            subset=self.subset, inplace=True,
            rmsd_cutoff=1.0, max_cycles=100
        )
        npdb1.write( self.superposed_file )
    def _post_exec( self ):
        self._make_provi_file(
            pdb_file1=self.relpath( self.pdb_file1 ),
            pdb_file2=self.relpath( self.pdb_file2 ),
            superposed_file=self.relpath( self.superposed_file ),
        )





def numpdb_test( pdb_file ):
    with Timer("read/parse pdb plain"):
        numpdb.NumPdb( pdb_file, features={
            "phi_psi": False, 
            "sstruc": False, 
            "backbone_only": True,
            "detect_incomplete": False
        })
    with Timer("read/parse pdb incomplete"):
        numpdb.NumPdb( pdb_file, features={"phi_psi": False, "sstruc": False, "backbone_only": False, "detect_incomplete": False} )
    with Timer("read/parse pdb phi/psi"):
        numpdb.NumPdb( pdb_file, features={"sstruc": False} )
    with Timer("read/parse pdb"):
        npdb = numpdb.NumPdb( pdb_file )
    # with Timer("resno iter2"):
    #     for numa_list in npdb.iter_resno2( 6, chain="A", resno=[1,22] ):
    #         print [ numa["resno"][0] for numa in numa_list ]
    # with Timer("dist"):
    #     print npdb.dist( {"chain":"A"}, {"chain":"B"} )
    with Timer("access phi/psi"):
        print np.nansum( npdb['phi'] ), np.nansum( npdb['psi'] )
    with Timer("access phi/psi, altloc"):
        print np.nansum( npdb.get('phi', altloc=[" ", "A"] ) )
    with Timer("sstruc iter"):
        for numa in npdb.iter_sstruc():
            pass
    with Timer("sequence"):
        print npdb.sequence(chain="A", resno=[1,100])
    
    npdb = numpdb.NumPdb( pdb_file, features={
        "phi_psi": False, 
        "sstruc": False, 
        "backbone_only": False,
        "detect_incomplete": False
    })
    with Timer("plain resno iter"):
        for numa in npdb._iter_resno():
            pass
    with Timer("write pdb"):
        npdb.write( "test.pdb" )

    npdb = numpdb.NumPdb( pdb_file, features={
        "phi_psi": False, 
        "sstruc": False, 
        "backbone_only": False,
        "detect_incomplete": True
    })
    with Timer("resno iter new"):
        for numa in npdb.iter_resno():
            pass


class NumpdbTest( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    def func( self ):
        numpdb_test( self.pdb_file )

