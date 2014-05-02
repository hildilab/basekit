import os
import unittest

from basekit.align import Muscle

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "align", *dir_name )


class MuscleTestCase( unittest.TestCase ):
    def setUp( self ):
        self.muscle = Muscle(
            data( "test.fasta" ),
            output_dir=tmp( "single_test" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.muscle()
        self.assertEquals( self.muscle.check( full=True ), "Ok" )

