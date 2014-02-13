import os
import unittest
import collections

from basekit.utils.tool import _, CmdTool, PyTool, RecordsMixin

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "utils_tool", *dir_name )



class SimpleCmdTool( CmdTool ):
    cmd = [ "echo", "Hello world!" ]

class CmdToolTestCase( unittest.TestCase ):
    def setUp( self ):
        self.tool = SimpleCmdTool( 
            output_dir=tmp( "simple_test" ),
            run=False
        )
    def test_check( self ):
        self.tool()
        self.assertEquals( self.tool.check( full=True ), "Ok" )
    def test_log( self ):
        self.tool()
        with open( self.tool.stdout_file, "r" ) as fp:
            log = fp.read()
        self.assertEquals( log, "Hello world!\n" )



SimpleRecord = collections.namedtuple("SimpleRecord", [ "name" ])

class SimpleRecordsTool( PyTool, RecordsMixin ):
    args = [
        _( "myid", type="text" )
    ]
    RecordsClass = SimpleRecord
    def _init( self, *args, **kwargs ):
        self._init_records( self.myid )
    def func( self ):
        pass

class RecordsToolTestCase( unittest.TestCase ):
    def setUp( self ):
        self.tool = SimpleRecordsTool(
            "foobar",
            output_dir=tmp( "records_test" ),
            run=False
        )
    def test_check( self ):
        self.tool()
        self.assertEquals( self.tool.check( full=True ), "Ok" )
