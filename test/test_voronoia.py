import os
import unittest

from basekit.voronoia import Voronoia

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )


# cd ./test/
# python -m unittest discover


class VoronoiaTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            os.path.join( DATA_DIR, "1T7H_A.pdb" ),
            ex=0.2,
            output_dir=os.path.join( TMP_DIR, "voronoia" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.voro()
        self.assertEquals( self.voro.check( full=True ), "Ok" )



class VoronoiaParallelTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            " ".join([
                os.path.join( DATA_DIR, "1T7H_A.pdb" ),
                os.path.join( DATA_DIR, "2WFU_A.pdb" ),
            ]),
            parallel="list",
            ex=0.2,
            output_dir=os.path.join( TMP_DIR, "voronoia" ),
            run=False 
        )
    def test_check( self ):
        self.voro()
        for t in self.voro.tool_results:
            self.assertEquals( t.check( full=True ), "Ok" )


