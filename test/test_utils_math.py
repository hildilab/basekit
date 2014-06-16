import unittest
import collections, os

import numpy as np
import basekit.utils.numpdb as numpdb
from basekit.utils.math import hclust, HCLUST_PACKAGE, rmatrixu, tmscore

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

class HclustTestCase( unittest.TestCase ):
    def test_hclust( self ):
        data = [
            [ 1, 1 ], [ 1, 2 ], [ 2, 1 ],
            [ 10, 10 ], [ 10, 11 ], [ 11, 10 ], [ 11, 11 ],
            [ 100, 101 ], [ 101, 100 ],
        ]
        clust = hclust( data, 2 )
        self.assertSequenceEqual(
            sorted( clust.values() ),
            sorted( collections.defaultdict( list, [
                ( 1, [[100, 101], [101, 100]] ),
                ( 2, [[1, 1], [1, 2], [2, 1]] ),
                ( 3, [[10, 10], [10, 11], [11, 10], [11, 11]] )
            ]).values() )
        )
    def test_hclust_package( self ):
        self.assertIn( HCLUST_PACKAGE, [ "fastcluster", "scipy" ] )
        self.assertEqual( 
            HCLUST_PACKAGE,
            "fastcluster", 
            "fastcluster not imported, using scipy.cluster.hierarchy instead" 
        )
        

class RotateTestCase( unittest.TestCase ):
    def test_rmatrixu( self ):
        positions = np.array([[ -1, 2, 0 ]])
        v = np.array([ 0, 2, 0 ])
        rotmat = rmatrixu( v, np.deg2rad( 180 ) )
        positions = ( np.dot( rotmat, positions.T ) ).T
        np.testing.assert_almost_equal( 
            positions, np.array([[ 1, 2, 0 ]]) 
        )

class TmscoreTestCase( unittest.TestCase ):
    def test_tmscore( self ):
        npdb1=numpdb.NumPdb( data("trajec_1.pdb") )
        coords1=npdb1.get('xyz', atomname='CA')
        npdb2=numpdb.NumPdb( data("trajec_2.pdb") )
        coords2=npdb2.get('xyz', atomname='CA')
        
        result=tmscore(coords1, coords2)
        np.testing.assert_almost_equal( 
            result, 0.780818185612
        )
