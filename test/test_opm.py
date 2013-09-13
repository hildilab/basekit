import os
import unittest

from basekit.opm import PPM_URL, _parse_ppm, parse_planes

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

# cd ./test/
# python -m unittest discover


class PpmParseTestCase( unittest.TestCase ):
    def setUp( self ):
        with open( data( "1a11_ppm.html" ), "r" ) as fp:
            self.html = fp.read()
    def test_parse( self ):
        pdb_url, info_dict = _parse_ppm( self.html )
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



