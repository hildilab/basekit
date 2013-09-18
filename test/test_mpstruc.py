import os
import unittest
import collections

from basekit.mpstruc import MpstrucDb

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )


# cd ./test/
# python -m unittest discover


class MpstrucDbTestCase( unittest.TestCase ):
    def setUp( self ):
        self.mpstruc = MpstrucDb()
    def test_list( self ):
        lst = self.mpstruc.list()
        dup = [ x for x, y in collections.Counter(lst).items() if y > 1 ]
        self.assertEquals( len( lst ) - len( dup ), 949 )
        self.assertListEqual( dup, [ '1BMF', '1A91' ] )
    def test_find_related_member( self ):
        info = self.mpstruc.info( "1AIG" )
        self.assertEquals( info['name'], "Photosynthetic Reaction Center" )
    def test_find_related_protein( self ):
        info = self.mpstruc.info( "2VQK" )
        self.assertEquals( info['name'], "Porin B monomer" )
    def test_find_member( self ):
        info = self.mpstruc.info( "1U19" )
        self.assertEquals( info['name'], "Rhodopsin" )
    def test_find_protein( self ):
        info = self.mpstruc.info( "1VGO" )
        self.assertEquals( info['name'], "Archaerhodopsin-2 (aR-2)" )


