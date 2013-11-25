import os
import json
import unittest

from basekit import utils, pdb
from basekit.utils import numpdb
from pdb import * 

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )
TMP_DIR2 = os.path.join( PARENT_DIR, "data","cionize" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "pdb", *dir_name )

def tmp2( *dir_name ):
    return os.path.join( TMP_DIR2, "cionize", *dir_name )


class RotamereTestCase( unittest.TestCase ):
    def test_make_rotamere( self ):
        npdb = numpdb.NumPdb( data( "testprot.pdb" ) )
        sele={ "resno": 20, "chain": 'A', "resname": 'TYR' }
        rotamere = pdb.make_rotamere( npdb, sele, 0 )
        self.assertAlmostEqual(
            rotamere.get( 'xyz', **sele )[-1][0],
            76.1959444
        )
        
class JoinSplittedTestCase( unittest.TestCase ):
    def test_join_splitted( self ):
        splitted = pdb.join_splitted( ['4GD1', '4GD2', '3R8S', '3R8T'], '', clear=True )
        self.assertEqual(
            splitted[0:20],
            'ATOM      1  P     A'
        )
        
class CionizeTestCase( unittest.TestCase ):
    def test_cionize( self ):
        ionized = pdb.Cionize(
            data( "3SN6.pdb" ),
            configfile=tmp2( "test.cfg" ),
            run=False,
            verbose=False)
        print ionized


class ClashTestCase ( unittest.TestCase ):
    def test_clash ( self ):
        npdb = numpdb.NumPdb("bestrotamerschain_n.pdb")
        hurz=pdb.find_all_clashes(npdb)

