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
@unittest.skipUnless(
    os.environ.get("LINKIT_DIR", False), 'LINDIT_DIR not set' )
class p2y12TestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "ssfe_p2y12" ), True )
        self.link_it = LinkIt(
            data( "ssfe/TMH_templates1_topmodel_rechain_p2y12.pdb" ),
            "88", "94", "TGPLR",
            output_dir=tmp( "ssfe_p2y12" ),
            run=False,
            verbose=True
        )

    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.link_it()
        self.assertEquals( self.link_it.check( full=True ), "Ok" )


@unittest.skipUnless(
    utils.path.which( 'wine' ), 'wine cmd not found' )
@unittest.skipUnless(
    os.environ.get("LINKIT_DIR", False), 'LINDIT_DIR not set' )
class p2y12NeuTestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "ssfe_p2y12_neu" ), True )
        self.multi_link_it = MultiLinkIt(
            data( "ssfe/TMH_templates1_topmodel_rechain_p2y12.pdb" ),
            [
                [ "53", "57", "RSK" ],
                [ "88", "94", "TGPLR" ],
                [ "128", "137", "PFKTSNPK" ],
                [ "161", "187", "LTNRQPRDKNVKKCSFLKSEFGLVW" ],
                [ "212", "230", "TKELYRSYVRTRGVGKV" ],
                [ "262", "276", "QTRDVFDCTAENT" ],
                [ "299", "302", "FL" ]
            ],
            names=[ "ICL1", "ECL1", "ICL2", "ECL2", "ICL3", "ECL3", "TM7-H8" ],
            output_dir=tmp( "ssfe_p2y12_neu" ),
            run=False,
            verbose=True
        )

    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.multi_link_it()
        self.assertEquals( self.multi_link_it.check( full=True ), "Ok" )


@unittest.skipUnless(
    utils.path.which( 'wine' ), 'wine cmd not found' )
@unittest.skipUnless(
    os.environ.get("LINKIT_DIR", False), 'LINDIT_DIR not set' )
class p2y12AltTestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "ssfe_p2y12_alt" ), True )
        self.multi_link_it = MultiLinkIt(
            data( "ssfe/TMH_templates1_topmodel_rechain_alt.pdb" ),
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
            output_dir=tmp( "ssfe_p2y12_alt" ),
            run=False,
            verbose=True
        )

    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.multi_link_it()
        self.assertEquals( self.multi_link_it.check( full=True ), "Ok" )


@unittest.skipUnless(
    utils.path.which( 'wine' ), 'wine cmd not found' )
@unittest.skipUnless(
    os.environ.get("LINKIT_DIR", False), 'LINDIT_DIR not set' )
class tshrAltTestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "ssfe_tshr_alt" ), True )
        self.multi_link_it = MultiLinkIt(
            data( "ssfe/TMH_templates1_topmodel_rechain_tshr_alt.pdb" ),
            [
                [ "443", "447", "YKL" ],
                [ "478", "491", "SEYYNHAIDWQT" ],
                [ "525", "534", "AMRLDRKI" ],
                [ "558", "578", "GISSYAKVSICLPMDTETP" ],
                [ "604", "618", "YITVRNPQYNPGD" ],
                [ "650", "656", "KPLIT" ],
                [ "679", "682", "IF" ]
            ],
            names=[ "ICL1", "ECL1", "ICL2", "ECL2", "ICL3", "ECL3", "TM7-H8" ],
            output_dir=tmp( "ssfe_tshr_alt" ),
            run=False,
            verbose=True
        )

    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.multi_link_it()
        self.assertEquals( self.multi_link_it.check( full=True ), "Ok" )


@unittest.skipUnless(
    utils.path.which( 'wine' ), 'wine cmd not found' )
@unittest.skipUnless(
    os.environ.get("LINKIT_DIR", False), 'LINDIT_DIR not set' )
class tshrNeuTestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "ssfe_tshr_neu" ), True )
        self.multi_link_it = MultiLinkIt(
            data( "ssfe/TMH_templates1_topmodel_rechain_tshr.pdb" ),
            [
                [ "443", "447", "YKL" ],
                [ "478", "491", "SEYYNHAIDWQT" ],
                [ "525", "534", "AMRLDRKI" ],
                [ "558", "578", "GISSYAKVSICLPMDTETP" ],
                [ "604", "618", "YITVRNPQYNPGD" ],
                [ "650", "656", "KPLIT" ],
                [ "679", "682", "IF" ]
            ],
            names=[ "ICL1", "ECL1", "ICL2", "ECL2", "ICL3", "ECL3", "TM7-H8" ],
            output_dir=tmp( "ssfe_tshr_neu" ),
            run=False,
            verbose=True
        )

    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.multi_link_it()
        self.assertEquals( self.multi_link_it.check( full=True ), "Ok" )


@unittest.skipUnless(
    utils.path.which( 'wine' ), 'wine cmd not found' )
@unittest.skipUnless(
    os.environ.get("LINKIT_DIR", False), 'LINDIT_DIR not set' )
class LinkitDensTestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "ribosomexample" ), True )
        self.linkitdens = LinkItDensity(
            data( "ribosomexample.pdb"),
            data(  "ribocut4a.mrc"),
            "154:C", "164:C", "EDKVEGYKK", 4.2,
            max_loops=100,
            output_dir=tmp( "ribosomexample" ),
            run=False,
            verbose=True
        )

    def test_check( self ):
        self.linkitdens()
        self.assertEquals( self.linkitdens.check( full=True ), "Ok" )


@unittest.skipUnless(
    utils.path.which( 'wine' ), 'wine cmd not found' )
@unittest.skipUnless(
    os.environ.get("LINKIT_DIR", False), 'LINDIT_DIR not set' )
class LinkItTestCase( unittest.TestCase ):
    def setUp( self ):
        shutil.rmtree( tmp( "3dqb_ICL3" ), True )
        self.link_it = LinkIt(
            data( "3DQB.pdb" ),
            "231:A", "248:A", "EAAAQQQESATTQKAE",
            max_loops=100,
            output_dir=tmp( "3dqb_ICL3" ),
            run=False,
            verbose=True
        )

    @unittest.skipIf( '-quick' in sys.argv, 'Long running' )
    def test_check( self ):
        self.link_it()
        self.assertEquals( self.link_it.check( full=True ), "Ok" )
