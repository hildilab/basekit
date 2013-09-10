import os
import unittest

from test import TMP_DIR, DATA_DIR

from basekit.voronoia import Voronoia


class VoronoiaTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            os.path.join( DATA_DIR, "1crn.pdb" ),
            ex=0.2,
            output_dir=os.path.join( TMP_DIR, "voronoia" ),
            run=False 
        )
    def test_check( self ):
        self.voro()
        self.assertEquals( self.voro.check( full=True ), "Ok" )



class VoronoiaParallelTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            " ".join([
                os.path.join( DATA_DIR, "1crn.pdb" ),
                os.path.join( DATA_DIR, "1u19.pdb" ),
            ]),
            parallel="list",
            ex=0.2,
            output_dir=os.path.join( TMP_DIR, "voronoia" ),
            run=False 
        )
    def test_check( self ):
        self.voro()
        for t in self.voro.results:
            self.assertEquals( t.check( full=True ), "Ok" )


