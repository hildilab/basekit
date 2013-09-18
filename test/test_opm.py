import os
import unittest

from basekit import utils
from basekit.opm import PPM_URL, _parse_opm_info, _parse_ppm, parse_planes

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "opm", *dir_name )


class OpmInfoParseTestCase( unittest.TestCase ):
    def test_parse( self ):
        with open( data( "1U19_opm_info.html" ), "r" ) as fp:
            info_dict = _parse_opm_info( fp.read() )
        self.assertDictEqual(
            info_dict,
            {
                'class': 'Alpha-helical polytopic',
                'delta_g': -71.9,
                'family': 'G-protein coupled receptors, family A',
                'localization': 'Eukaryotic plasma membrane',
                'related_ids': ['1F88', '1HZX', '1L9H', '2G87'],
                'species': 'taurus',
                'superfamily': 'Rhodopsin-like receptors and pumps',
                'type': 'Transmembrane'
            }
        )
    def test_parse_related( self ):
        with open( data( "3DQB_opm_info.html" ), "r" ) as fp:
            info_dict = _parse_opm_info( fp.read() )
        self.assertDictEqual(
            info_dict,
            {
                "representative": "3PQR"
            }
        )


class PpmParseTestCase( unittest.TestCase ):
    def tearDown( self ):
        utils.path.remove( "ppm_error.txt" )
    def test_parse( self ):
        with open( data( "1a11_ppm.html" ), "r" ) as fp:
            pdb_url, info_dict = _parse_ppm( fp.read() )
        self.assertEqual(
            pdb_url,
            PPM_URL + "pdb_upload/1a11_asmout.pdb"
        )
        self.assertDictEqual(
            info_dict,
            {
                "delta_g": -17.8
            }
        )
    def test_error( self ):
        with self.assertRaises( Exception ) as context:
            with open( data( "ppm_error.html" ), "r" ) as fp:
                pdb_url, info_dict = _parse_ppm( fp.read() )
        self.assertEqual(
            context.exception.message,
            "./transversion < ./tmp/1f442ced987451199522e5ae0c56c1eb.in > "
            "./tmp/1f442ced987451199522e5ae0c56c1eb.out"
        )
    def test_error2( self ):
        with self.assertRaises( Exception ) as context:
            with open( data( "ppm_error2.html" ), "r" ) as fp:
                pdb_url, info_dict = _parse_ppm( fp.read() )
        self.assertEqual(
            context.exception.message,
            "Too many residues"
        )


class ParsePlanesTestCase( unittest.TestCase ):
    def test_parse_hetatom_dum( self ):
        mplanes = parse_planes( data( "3DQB_opm.pdb" ) )
        self.assertListEqual(
            mplanes.tolist(),
            [
                [
                    [-26., -8., -15.95],
                    [-26., -6., -15.95],
                    [-24., -12., -15.95]
                ],
                [
                    [-26., -8., 15.95],
                    [-26., -6., 15.95],
                    [-24., -12., 15.95]
                ]
            ]
        )
    def test_parse_atom_dum( self ):
        mplanes = parse_planes( data( "4JR9_opm.pdb" ) )
        self.assertListEqual(
            mplanes.tolist(),
            [
                [
                    [-30.0, -4.0, -14.8], 
                    [-30.0, -2.0, -14.8], 
                    [-28.0, -12.0, -14.8]
                ],
                [
                    [-30.0, -4.0, 14.8],
                    [-30.0, -2.0, 14.8],
                    [-28.0, -12.0, 14.8]
                ]
            ]
        )
    def test_no_mplane_found_exception( self ):
        with self.assertRaises( Exception ) as context:
            parse_planes( data( "1U19.pdb" ) )
        self.assertEqual(
            context.exception.message,
            "could not find plane coordinates"
        )
    def test_parse_single_plane( self ):
        mplanes = parse_planes( data( "1BA4_opm.pdb" ) )
        self.assertListEqual(
            mplanes.tolist(),
            [
                [
                    [-22.0, -6.0, 15.4],
                    [-22.0, -4.0, 15.4],
                    [-20.0, -10.0, 15.4]
                ]
            ]
        )



