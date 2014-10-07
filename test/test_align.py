import os
import unittest

from basekit.align import Muscle, TheseusMakeFasta

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
        self.muscle2 = Muscle(
            data( "test.fasta" ),
            pdb_files=[
                data( "1U19.pdb" ),
                data( "3SN6.pdb" ),
            ],
            output_dir=tmp( "single_test" ),
            mapfile=True,
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.muscle()
        self.muscle2()
        self.assertEquals( self.muscle.check( full=True ), "Ok" )
        self.assertEquals( self.muscle2.check( full=True ), "Ok" )

class TheseusMakeFastaTestCase( unittest.TestCase ):
    def setUp( self ):
        self.tmf = TheseusMakeFasta(
            [
                data( "1U19.pdb" ),
                data( "3SN6.pdb" ),
            ],
            output_dir=tmp( "single_test" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.tmf()
        self.assertEquals( self.tmf.check( full=True ), "Ok" )