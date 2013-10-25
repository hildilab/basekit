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
            "configuration": False,
            "info": True
        })
        print "\n" + json.dumps( npdb._info, indent=4 )
        print "\n"
        for ss in npdb._sstruc:
            print ss
        print "\n"
        for ss in npdb.iter_sstruc():
            print ss._atoms["resno"].min(), ss._atoms["resno"].max()


class RotamereTestCase( unittest.TestCase ):
    def test_make_rotamere( self ):
        npdb = numpdb.NumPdb( data( "testprot.pdb" ))
        sele={ "resno": 21, "chain": "A", "resname": "VAL" }
        no = numpdb.get_rotno ( sele["resname"] )
        for i in range(0, 1):#no):
            print '###NEXT ROUND###', i
            rotamere = numpdb.make_rotamere( npdb, sele, i )
        

        #positions = np.array([[ -1, 2, 0 ]])
        #v = np.array([ 0, 1, 0 ])
        #rotmat = rmatrixu( v, np.deg2rad( 180 ) )
        #positions = ( np.dot( rotmat, positions.T ) ).T
        #np.testing.assert_almost_equal( 
        #    positions, np.array([[ 1, 2, 0 ]]) 
        #)
        
        
        
        