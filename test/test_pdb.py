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

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "pdb", *dir_name )


class RotamereTestCase( unittest.TestCase ):
    def test_make_rotamere( self ):
        npdb = numpdb.NumPdb( data( "testprot.pdb" ) )
        sele={ "resno": 20, "chain": 'A', "resname": 'TYR' }
        rotamere = pdb.make_rotamere( npdb, sele, 0 )
        self.assertAlmostEqual(
            rotamere.get( 'xyz', **sele )[-1][0],
            76.1959444
        )