import os
import unittest

from basekit.superpose import Theseus

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "superpose", *dir_name )


class SuperposeTestCase( unittest.TestCase ):
    def setUp( self ):
        self.theseus = Theseus(
            [data( "1U19.pdb" ), data( "1CRN.pdb" )],
            prefix="test",
            output_dir=tmp( "single_test" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.theseus()
        self.assertEquals( self.theseus.check( full=True ), "Ok" )
