import os
import unittest
import collections

from basekit import utils
from basekit.mpstruc import MpstrucDb

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )


# cd ./test/
# python -m unittest discover

@unittest.skipUnless(
    utils.path.which( 'MPSTRUC_ALPHAHELICAL_PATH' ), 'MPSTRUC_ALPHAHELICAL_PATH not set' )
class MpstrucDbTestCase( unittest.TestCase ):
    def setUp( self ):
        self.mpstruc = MpstrucDb()
    def test_list( self ):
        lst = self.mpstruc.list()
        dup = [ x for x, y in collections.Counter(lst).items() if y > 1 ]
        # FIXME the length changes when the MPstruc db changes
        # self.assertEquals( len( lst ) - len( dup ), 949 )
        self.assertListEqual( dup, [ '1BMF', '1A91' ] )
    def test_find_related_member( self ):
        info = self.mpstruc.info( "1AIG" )
        self.assertEquals( info['name'], "Photosynthetic Reaction Center" )
        self.assertDictEqual(
            info,
            {
                'group': None,
                'master': '4RCR',
                'name': 'Photosynthetic Reaction Center',
                'related': [ '1AIG' ],
                'species': 'Rhodobacter sphaeroides (dark state)',
                'subgroup': None
            }
        )
    def test_find_related_protein( self ):
        info = self.mpstruc.info( "2VQK" )
        self.assertItemsEqual( info[ "related" ], ['2VQH', '2VQK', '2VQL'] )
        del info[ "related" ]
        self.assertDictEqual(
            info,
            {
                'group': 'TRANSMEMBRANE PROTEINS: ALPHA-HELICAL',
                'master': None,
                'name': 'Porin B monomer',
                'species': 'Corynebacterium glutamicum',
                'subgroup': 'Outer Membrane Proteins'
            }
        )
    def test_find_member( self ):
        info = self.mpstruc.info( "1U19" )
        self.assertDictEqual(
            info,
            {
                'group': 'TRANSMEMBRANE PROTEINS: ALPHA-HELICAL',
                'master': '1F88',
                'name': 'Rhodopsin',
                'related': [],
                'species': 'Bovine Rod Outer Segment',
                'subgroup': 'G Protein-Coupled Receptors (GPCRs)'
            }
        )
    def test_find_protein( self ):
        info = self.mpstruc.info( "1VGO" )
        self.assertEquals( info['name'], "Archaerhodopsin-2 (aR-2)" )
        self.assertDictEqual(
            info,
            {
                'group': 'TRANSMEMBRANE PROTEINS: ALPHA-HELICAL',
                'master': None,
                'name': 'Archaerhodopsin-2 (aR-2)',
                'related': [],
                'species': 'Haloroubrum sp. aus-2',
                'subgroup': 'Bacterial and Algal Rhodopsins'
            }
        )
    def test_format( self ):
        info = self.mpstruc.info( "3SN6" )
        self.assertEquals(
            info['name'],
            "beta2 adrenergic receptor-Gs protein complex"
        )
    def test_format2( self ):
        info = self.mpstruc.info( "3DQB" )
        self.assertEquals(
            info['name'],
            "Rhodopsin, Ops*-GalphaCT peptide complex"
        )



