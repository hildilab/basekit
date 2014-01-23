import os
import sys
import shutil
import unittest

from basekit import utils
from basekit.linker import LinkIt, MultiLinkIt, LINKIT_CMD, LinkItDensity

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "linker", *dir_name )


@unittest.skipUnless( 
        utils.path.which( 'wine' ), 'wine cmd not found' )
class LinkItTestCase( unittest.TestCase ):
    def setUp( self ):
        self.link_it = LinkIt(
            data( "ssfe_p2y12.pdb" ),
            "88", "94", "TGPLR",
            output_dir=tmp( "ssfe_p2y12" ),
            run=False,
            verbose=False
        )
    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.link_it()
        self.assertEquals( self.link_it.check( full=True ), "Ok" )



@unittest.skipUnless( 
        utils.path.which( 'wine' ), 'wine cmd not found' )
class MultiLinkItTestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "ssfe_p2y12_multi" ), True )
        self.multi_link_it = MultiLinkIt(
            data( "ssfe_p2y12.pdb" ),
            [
                [ "53", "57", "RSK" ],
                [ "88", "94", "TGPLR" ],
                [ "128", "137", "PFKTSNPK" ],
                [ "161", "186", "LTNRQPRDKNVKKCSFLKSEFGLV" ],
                [ "212", "230", "TKELYRSYVRTRGVGKV" ],
                [ "262", "276", "QTRDVFDCTAENT" ],
                [ "299", "302", "FL" ]
            ],
            names=[ "ICL1", "ECL1", "ICL2", "ECL2", "ICL3", "ECL3", "TM7-H8" ],
            output_dir=tmp( "ssfe_p2y12_multi" ),
            run=False,
            verbose=False
        )
    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.multi_link_it()
        self.assertEquals( self.multi_link_it.check( full=True ), "Ok" )


        
@unittest.skipUnless( 
        utils.path.which( 'wine' ), 'wine cmd not found' )
class LinkitDensTestCase( unittest.TestCase ):
    def setUp( self ):
        self.linkitdens= LinkItDensity(
            data( "ribosomexample.pdb", "ribocut4a.mrc"),
            "154:C", "164:C", "EDKVEGYKK", "4.2",
            output_dir=tmp( "ribosomexample" ),
            run=False,
            verbose=False
        )
    def set_check( self ):
        self.linkitdens()
        self.assertEquals( self.linkitdens( full=True ), "Ok" )

