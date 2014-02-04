import unittest
import os

from basekit.utils.mrc import get_mrc_header

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )


class MrcHeaderTestCase( unittest.TestCase ):
    def test_header( self ):
        mrc_h = get_mrc_header( data( "ribocut4a.mrc" ) )
        self.assertEqual( mrc_h.cmap, 'MAP ' )
        self.assertEqual( mrc_h.stamp, 16708 )
        self.assertEqual( mrc_h.xorg, 0 )
        self.assertEqual( mrc_h.yorg, 0 )
        self.assertEqual( mrc_h.zorg, 0 )
        self.assertEqual( mrc_h.amax, 35971.1953125 )
        self.assertEqual( mrc_h.amin, -19664.234375 )
