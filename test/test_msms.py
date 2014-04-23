import os
import unittest

from basekit.msms import Msms

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "msms", *dir_name )


class MsmsTestCase( unittest.TestCase ):
    def setUp( self ):
        self.msms = Msms(
            data( "1CRN.pdb" ),
            all_components=True,
            output_dir=tmp( "single_test" ),
            verbose=False
        )
    def test_check( self ):
        self.msms()
        self.assertEquals( len(self.msms.obj_list), 1 )
        self.assertEquals( self.msms.check( full=True ), "Ok" )

