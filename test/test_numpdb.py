import os
import json
import unittest

from basekit import utils
from basekit.utils import numpdb

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "numpdb", *dir_name )


class NumpdbHeaderParseTestCase( unittest.TestCase ):
    def test_secondary_structure( self ):
        npdb = numpdb.NumPdb( data( "1CRN.pdb" ), {
            "phi_psi": False,
            "sstruc": True,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False
        })
        self.assertEquals(
            npdb._sstruc[0].resname1, "THR"
        )
    def test_info( self ):
        npdb = numpdb.NumPdb( data( "1CRN.pdb" ), {
            "phi_psi": False,
            "sstruc": True,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False,
            "info": True
        })
        self.assertDictEqual(
            npdb._info,
            {
                "title": "WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC "
                    "RESOLUTION. PENTAGON RINGS OF WATER MOLECULES IN CRYSTALS "
                    "OF CRAMBIN", 
                "obsolete": [], 
                "model_type": {}, 
                "experiment": "X-RAY DIFFRACTION", 
                "keywords": [
                    "PLANT SEED PROTEIN", 
                    "PLANT PROTEIN"
                ], 
                "splited_entry": [], 
                "resolution": 1.5
            }
        )
    def test_missing( self ):
        npdb = numpdb.NumPdb( data( "3SN6.pdb" ), {
            "phi_psi": False,
            "sstruc": True,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False,
            "detect_missing": True
        })
        self.assertEqual(
        str(npdb._missing[7]),
        "MissingRecord(type='Residue', model='    ', "+
        "resname='MET', chain=' A', ssseqi=1, identifier='  ', "+
        "atoms=None)")



class NumseleTestCase( unittest.TestCase ):
    def test_chain( self ):
        self.assertDictEqual(
            numpdb.numsele( ":A" ),
            { "chain": "A", "resno": None, "atomname": None }
        )
        self.assertDictEqual(
            numpdb.numsele( ": " ),
            { "chain": " ", "resno": None, "atomname": None }
        )
    def test_chain_too_long( self ):
        with self.assertRaises( Exception ) as context:
            numpdb.numsele( ":AX" )
        self.assertEqual(
            context.exception.message,
            "chain identifier must be one character"
        )
        with self.assertRaises( Exception ) as context:
            numpdb.numsele( ":AX.CA" )
        self.assertEqual(
            context.exception.message,
            "chain identifier must be one character"
        )
    def test_resno( self ):
        self.assertDictEqual(
            numpdb.numsele( "312" ),
            { "chain": None, "resno": 312, "atomname": None }
        )
    def test_resno_range( self ):
        self.assertDictEqual(
            numpdb.numsele( "110-311" ),
            { "chain": None, "resno": [ 110, 311 ], "atomname": None }
        )
        self.assertDictEqual(
            numpdb.numsele( "110-312:C" ),
            { "chain": "C", "resno": [ 110, 312 ], "atomname": None }
        )
        self.assertDictEqual(
            numpdb.numsele( "110-313.N" ),
            { "chain": None, "resno": [ 110, 313 ], "atomname": "N" }
        )
    def test_atomname( self ):
        self.assertDictEqual(
            numpdb.numsele( ".CA" ),
            { "chain": None, "resno": None, "atomname": "CA" }
        )
    def test_atomname_too_long( self ):
        with self.assertRaises( Exception ) as context:
            numpdb.numsele( ":A.foobar" )
        self.assertEqual(
            context.exception.message,
            "atomname must be one to four characters"
        )
    def test_combination( self ):
        self.assertDictEqual(
            numpdb.numsele( "1:B.CA" ),
            { "chain": "B", "resno": 1, "atomname": "CA" }
        )
        self.assertDictEqual(
            numpdb.numsele( "1: .CA" ),
            { "chain": " ", "resno": 1, "atomname": "CA" }
        )
    
