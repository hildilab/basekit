import os
import unittest


try:
    from basekit import moderna_tools

    DIR = os.path.split( os.path.abspath( __file__ ) )[0]
    PARENT_DIR = os.path.split( DIR )[0]
    DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )

    def data( file_name ):
        return os.path.join( DATA_DIR, file_name )


    class ModernaTestCase( unittest.TestCase ):
        def test_format( self ):
            load = moderna_tools.load_templates( data( "3SN6.pdb" ) )
            self.assertEquals(
                str(load),
                "<RNA structure; chain 'A'; 349 residues>"
            )

except ImportError, e:
    class TestMissingDependency(unittest.TestCase):
        @unittest.skip( e )
        def test_fail():
            pass



