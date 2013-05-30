#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division



import sys
import os
import argparse
import functools
import itertools
import inspect
import json
import string
import collections

from basekit.utils import try_int, get_index, boolean, working_directory
from basekit.utils.timer import Timer
from basekit.utils.job import run_command


TIMEOUT_CMD = "timeout"



def get_type( params ):
    if params.get('fixed'):
        return float
    elif params["type"]=="slider":
        return int
    elif params["type"] in [ "file", "text", "select" ]:
        return str
    elif params["type"]=="checkbox":
        return boolean
    else:
        return str

def make_parser( Tool, parser=None ):
    if not parser:
        parser = argparse.ArgumentParser( description=Tool.__doc__, epilog="basekit", add_help=False )
    for name, params in Tool.args.iteritems():
        option = '--%s'%name if "default_value" in params else name
        default = params.get( "default_value", None )
        type = get_type( params )
        parser.add_argument( option, type=type, default=default, help=params.get("help") )
    if not Tool.no_output:
        group = parser.add_argument_group( title="general arguments" )
        group.add_argument( '-o', metavar='OUTPUT_DIR', type=str, default="./" )
        group.add_argument( '-t', metavar='TIMEOUT', type=int, default=0 )
        group.add_argument( '-v', '--verbose', action='store_true' )
        group.add_argument( '-c', '--check', action='store_true' )
        group.add_argument( '-h', '--help', action="help", help="show this help message and exit" )
    return parser

def parse_args( Tool, kwargs=None ):
    if not kwargs:
        parser = make_parser( Tool )
        kwargs = vars( parser.parse_args() )
    args = []
    for name, params in Tool.args.iteritems():
        if "default_value" not in params:
            args.append( kwargs.pop( name ) )
    return args, kwargs

def parse_subargs( tools, description=None ):
    if description:
        description += " The tools are accessible by the subcommands given below."
    else:
        description = "A collection of tools, accessible by the subcommands given below."
    parser = argparse.ArgumentParser( description=description )
    subparsers = parser.add_subparsers( title='subcommands' )
    for name, Tool in tools.iteritems():
        subp = subparsers.add_parser( name, description=Tool.__doc__, epilog="basekit", add_help=False )
        make_parser( Tool, subp )
        subp.set_defaults(Tool=Tool)
    pargs = vars( parser.parse_args() )
    Tool = pargs.pop("Tool")
    args, kwargs = parse_args( Tool, kwargs=pargs )
    return Tool, args, kwargs



class ToolMetaclass(type):
    def __init__(cls, name, bases, dct):
        if not "no_output" in dct:
            cls.no_output = False

        args = collections.OrderedDict()
        for a in dct.get( "args", [] ):
            args[ a.pop("name") ] = a
        cls.args = args



class Tool( object ):
    __metaclass__ = ToolMetaclass
    def __init__( self, *args, **kwargs ):
        self.name = self.__class__.__name__.lower()

        self.timeout = kwargs.get("timeout", None)
        self.fileargs = kwargs.get("fileargs", False)
        self.verbose = kwargs.get("verbose", False)
        self.output_dir = os.path.abspath( kwargs.get("output_dir", ".") ) + os.sep

        if not self.no_output:
            if not os.path.exists( self.output_dir ):
                os.makedirs( self.output_dir )

            self.args_file = os.path.join( self.output_dir, "%s.json" % self.name )
            if self.fileargs:
                with open( self.args_file, "r" ) as fp:
                    args, kwargs = json.load( fp )
            else:
                with open( self.args_file, "w" ) as fp:
                    json.dump( ( args, kwargs ), fp, indent=4 )
        
        self._init( *args, **kwargs )

        if kwargs.get("run", True) and not kwargs.get("check", False):
            self.__run()
    def __run( self ):
        with working_directory( self.output_dir ):
            self._pre_exec()
            self._run()
            self._post_exec()
    def __call__( self ):
        self.__run()
        return self
    def _run( self ):
        pass
    def _pre_exec( self ):
        pass
    def _post_exec( self ):
        pass
    def __check_file( self, f ):
        f = os.path.abspath(f)
        try:
            return os.path.isfile(f) and os.path.getsize(f)>0
        except:
            return False
    def check( self, full=False ):
        if hasattr( self, "output_files" ):
            with working_directory( self.output_dir ):
                isfile = map( self.__check_file, self.output_files )
                if full:
                    return list( isfile )
                else:
                    return all( isfile )
        else:
            return True
    def __str__( self ):
        status = "ok" if self.check() else "failed"
        return "%s status: %s" % ( self.name, status )
    def relpath( self, path, no_ext=False ):
        if no_ext:
            path = os.path.splitext( path )[0]
        return os.path.relpath( path, self.output_dir )
    def outpath( self, file_name ):
        return os.path.join( self.output_dir, file_name )
    def abspath( self, path ):
        return os.path.abspath( path )
    def subdir( self, directory ):
        return os.path.join( self.output_dir, directory )




class PyTool( Tool ):
    def __init__( self, *args, **kwargs ):
        super(PyTool, self).__init__( *args, **kwargs )
        if not hasattr( self, "func" ):
            raise Exception("A PyTool needs a 'func' attribute")
    def _run( self ):
        return self.func()


class CmdTool( Tool ):
    def __init__( self, *args, **kwargs ):
        super(CmdTool, self).__init__( *args, **kwargs )
        if not hasattr( self, "cmd" ):
            raise Exception("A CmdTool needs a 'cmd' attribute")
    def _run( self ):
        cmd = self.cmd
        if self.timeout:
            cmd = [ TIMEOUT_CMD, self.timeout ] + cmd
        log_file = "%s_cmd.log" % self.name
        ret = run_command( cmd, log=log_file, verbose=self.verbose )
        return ret




class ScriptMixin( object ):
    def _make_script_file( self, **values_dict ):
        tmpl_file = os.path.join( self.tmpl_dir, self.tmpl_file )
        with open( tmpl_file, "r" ) as fp:
            tmpl_str = fp.read()
        script_file = os.path.join( self.output_dir, self.tmpl_file )
        with open( script_file, "w" ) as fp:
            fp.write( string.Template( tmpl_str ).substitute( **values_dict ) )
        return script_file




def do_parallel( tool, files, args=None, kwargs=None, nworkers=None, run=True ):
    # !important - allows one to abort via CTRL-C
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    multiprocessing.log_to_stderr( logging.ERROR )
    
    if not nworkers: nworkers = multiprocessing.cpu_count()
    p = multiprocessing.Pool( nworkers )

    if not kwargs: kwargs = {}
    kwargs["run"] = run
    if not args: args = []
    data = p.map( functools.partial( tool, *args, **kwargs ), files ) # does partial work with a functor?
    p.close()
    p.join()

    return data






