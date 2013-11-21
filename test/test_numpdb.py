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
        npdb = numpdb.NumPdb( data( "3SN6.pdb" ), {#1CRN.pdb" ), {
            "phi_psi": False,
            "sstruc": True,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False,
            "info": True,
            "detect_missing": True
        })
        print "\n" + json.dumps( npdb._info, indent=4 )
        print "\n"
        for ss in npdb._sstruc:
            print ss
        print "\n"
        for ss in npdb.iter_sstruc():
            print ss._atoms["resno"].min(), ss._atoms["resno"].max()
        print "\n"
        for ss in npdb._missing:
            print ss
