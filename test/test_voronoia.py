import os
import unittest

from basekit.voronoia import Voronoia, make_ref, parse_vol

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "voronoia", *dir_name )


class VoronoiaTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            data( "1T7H_A.pdb" ),
            ex=0.3,
            output_dir=tmp( "single_test" ),
            run=False,
            verbose=False
        )
    def test_check( self ):
        self.voro()
        self.assertEquals( self.voro.check( full=True ), "Ok" )



class VoronoiaParallelTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            " ".join([ data( "1T7H_A.pdb" ), data( "2WFU_A.pdb" ) ]),
            parallel="list",
            ex=0.3,
            output_dir=tmp( "parallel_test" ),
            run=False 
        )
    def test_check( self ):
        self.voro()
        for t in self.voro.tool_results:
            self.assertEquals( t.check( full=True ), "Ok" )



class VoronoiaMakeRefTestCase( unittest.TestCase ):
    def setUp( self ):
        self.voro = Voronoia(
            " ".join([
                os.path.join( DATA_DIR, "1T7H_A.pdb" ),
                os.path.join( DATA_DIR, "2WFU_A.pdb" ),
            ]),
            parallel="list",
            ex=0.3,
            output_dir=os.path.join( TMP_DIR, "voronoia" ),
            make_reference=True,
            run=False 
        )
        self.voro()
        self.mr = make_ref( self.voro.tool_results )
    def test_check( self ):
        dic, dic2, log, pdb = self.mr
        self.assertNotEqual(dic, [])
        self.assertNotEqual(dic2, [])
        #self.assertNotEqual(log, '')

class VoronoiaParseVolTestCase( unittest.TestCase ):
    def setUp( self ):
        self.pv = parse_vol(
            os.path.join( DATA_DIR, "1T7H_A.vol" ),
            os.path.join( DATA_DIR, "1T7H_A.pdb" )
        )
    def test_check( self ):
        self.assertEqual( sum( self.pv['nrholes'].values() ), 0 )
        self.assertEqual(self.pv['holes'], [])
        self.assertEqual(sum( self.pv['packdens'].values()), 73.99856896262003)
        self.assertEqual(sum( self.pv['buried'].values()), 26)
        self.assertNotEqual(self.pv['zscore'], {})
        self.assertNotEqual(self.pv['protors'], {})
        #self.assertNotEqual(self.pv['log_list'], '')

