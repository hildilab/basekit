import unittest
import collections

from basekit.utils.math import hclust, HCLUST_PACKAGE


class HclustTestCase( unittest.TestCase ):
    def test_hclust( self ):
        data = [
            [ 1, 1 ], [ 1, 2 ], [ 2, 1 ],
            [ 10, 10 ], [ 10, 11 ], [ 11, 10 ],
            [ 100, 101 ], [ 101, 100 ],
        ]
        clust = hclust( data, 2 )
        self.assertEqual(
            clust,
            collections.defaultdict( list, [
                ( 1, [[100, 101], [101, 100]] ),
                ( 2, [[1, 1], [1, 2], [2, 1]] ),
                ( 3, [[10, 10], [10, 11], [11, 10]] )
            ])
        )
    def test_hclust_package( self ):
        self.assertIn( HCLUST_PACKAGE, [ "fastcluster", "scipy" ] )
        self.assertEqual( 
            HCLUST_PACKAGE,
            "fastcluster", 
            "fastcluster not imported, using scipy.cluster.hierarchy instead" 
        )
        


