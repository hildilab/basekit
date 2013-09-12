import os
import unittest

from basekit.opm import PPM_URL, _parse_ppm

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )


# cd ./test/
# python -m unittest discover


class PpmParseTestCase( unittest.TestCase ):
    def setUp( self ):
    	with open( os.path.join( DATA_DIR, "1a11_ppm.html" ), "r" ) as fp:
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

