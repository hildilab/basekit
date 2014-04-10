import unittest
import collections
import numpy as np


from basekit.utils.listing import merge_dic_list



class MergeTestCase( unittest.TestCase ):
    def test_merge_dic_list( self ):
        dic1={
            "a":{
                1:[1, 1.1, 1.11],
                2:[2,2.2],
                3:[2,2.2],
            }, "b":{
                1:{
                    "one":["o", "n", "e"],
                    "two":["t", "w", "o"]
                }
            }
        }
        dic2={
            "a":{
                1:[1, 1.1, 1.11, 1.111],
                2:{
                    "a":[1],
                    "b":[1,2]
                },
                4:{
                    "a":[1],
                    "b":[1,2]
                },
                5:[5] 
            }, "b":{
                1:{
                    "one":["o", "n", "e", "s"],
                    "two":["t"]
                }
            }, "c":{
                1:{
                    "I":["I", "II"],
                    "II":["I", "II"]
                }
            }
        }
        dic3={
            'a': {
                1: [1, 1.1, 1.11, 1, 1.1, 1.11, 1.111],
                2: {
                    'a': [1],
                    'b': [1, 2]
                },
                3: [2, 2.2],
                4: {
                    'a': [1],
                    'b': [1, 2]
                },
                5: [5]
            }, 'c': {
                1: {
                    'I': ['I', 'II'],
                    'II': ['I', 'II']
                }
            }, 'b': {
                1: {
                    'two': ['t', 'w', 'o', 't'],
                    'one': ['o', 'n', 'e', 'o', 'n', 'e', 's']
                }
            }
        }
        result=merge_dic_list(dic1, dic2)
        np.testing.assert_equal( 
            dic3, result 
        )
        
