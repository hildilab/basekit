import os
import unittest

from basekit.voronoia import Voronoia

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "voronoia", *dir_name )


class VoronoiaTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            data( "1T7H_A.pdb" ),
            ex=0.3,
            output_dir=tmp( "single_test" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.voro()
        self.assertEquals( self.voro.check( full=True ), "Ok" )



class VoronoiaParallelTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            " ".join([ data( "1T7H_A.pdb" ), data( "2WFU_A.pdb" ) ]),
            parallel="list",
            ex=0.3,
            output_dir=tmp( "parallel_test" ),
            run=False 
        )
    def test_check( self ):
        self.voro()
        for t in self.voro.tool_results:
            self.assertEquals( t.check( full=True ), "Ok" )


