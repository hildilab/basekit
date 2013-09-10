import unittest

from test import TMP_DIR

from basekit.utils.tool import CmdTool





class HelloWorldTool( CmdTool ):
    cmd = [ "echo", "-e", "Hello world!" ]


class CmdToolTestCase( unittest.TestCase ):
    def setUp( self ):
        self.tool = HelloWorldTool( output_dir=TMP_DIR, run=False )
    def test_check( self ):
        self.tool()
        self.assertEquals( self.tool.check( full=True ), "Ok" )
    def test_log( self ):
        self.tool()
        with open( self.tool.stdout_file, "r" ) as fp:
            log = fp.read()
        self.assertEquals( log, "Hello world!\n" )

