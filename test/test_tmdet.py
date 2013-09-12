import os
import unittest

from basekit.tmdet import PdbtmDb

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )


# cd ./test/
# python -m unittest discover


class PdbtmDbTestCase( unittest.TestCase ):
    def setUp( self ):
        self.pdbtm = PdbtmDb( os.path.join( DATA_DIR, "pdbtm_test.xml" ) )
    def test_find( self ):
        ret = self.pdbtm.find( "3a0b" )
        self.assertEquals( ret.attrib["ID"].lower(), "3a0b" )
    def test_list( self ):
        lst = self.pdbtm.list()
        self.assertEquals( lst,  ['2a06', '2a0d', '3a0b', '3a0h'] )

