import os
import unittest

from basekit.mapman import Mapman, BrixToMap

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "mapman", *dir_name )


class MapmanTestCase( unittest.TestCase ):
    def setUp( self ):
        self.mapman = Mapman(
            data( "4bs3.brix" ),
            "brix",
            output_dir=tmp( "single_test" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.mapman()
        self.assertEquals( self.mapman.check( full=True ), "Ok" )

class BrixToMapTestCase( unittest.TestCase ):
    def setUp( self ):
        self.brix2map = BrixToMap(
            data( "4bs3.brix" ),
            output_dir=tmp( "single_test" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.brix2map()
        self.assertEquals( self.brix2map.check( full=True ), "Ok" )
