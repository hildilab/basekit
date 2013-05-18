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
from string import Template
from collections import OrderedDict

from basekit.utils import try_int, get_index, boolean, working_directory
from basekit.utils.timer import Timer
from basekit.utils.job import run_command


TIMEOUT_CMD = "timeout"


def make_args( args ):
    _args = OrderedDict()
    for a in args:
        _args[ a.pop("name") ] = a
    return _args

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
        parser = argparse.ArgumentParser()
    for name, params in Tool.args.iteritems():
        option = '--%s'%name if "default_value" in params else name
        default = params.get( "default_value", None )
        type = get_type( params )
        parser.add_argument( option, type=type, default=default)    
    parser.add_argument( '-o', '--output_dir', type=str, default="./" )
    parser.add_argument( '-t', '--timeout', type=int, default=0 )
    parser.add_argument( '-v', '--verbose', action='store_true' )
    parser.add_argument( '-c', '--check', action='store_true' )
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

def parse_subargs( tools ):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers( title='subcommands' )
    for name, Tool in tools.iteritems():
        subp = subparsers.add_parser( name )
        make_parser( Tool, subp )
        subp.set_defaults(Tool=Tool)
    pargs = vars( parser.parse_args() )
    Tool = pargs.pop("Tool")
    args, kwargs = parse_args( Tool, kwargs=pargs )
    return Tool, args, kwargs





class Tool( object ):
    args = make_args([])
    def __init__( self, *args, **kwargs ):
        self.name = self.__class__.__name__.lower()

        self.timeout = kwargs.get("timeout", None)
        self.fileargs = kwargs.get("fileargs", False)
        self.verbose = kwargs.get("verbose", False)

        self.output_dir = os.path.abspath( kwargs.get("output_dir", ".") ) + os.sep
        if not os.path.exists( self.output_dir ):
            os.makedirs( self.output_dir )

        self.args_file = os.path.join( self.output_dir, "%s.json" % self.name )
        if self.fileargs:
            with open( self.args_file, "r" ) as fp:
                args, kwargs = json.load( fp )
        else:
            with open( self.args_file, "w" ) as fp:
                json.dump( ( args, kwargs ), fp, indent=4 )
        print args, kwargs
        self._init( *args, **kwargs )

        if kwargs.get("run", True) and not kwargs.get("check", True):
            self.__run()
    def __run( self ):
        with working_directory( self.output_dir ):
            self._pre_exec()
            self._run()
    def __call__( self ):
        self.__run()
        return self
    def _run( self ):
        pass
    def _pre_exec( self ):
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




class PyTool( Tool ):
    def __init__( self, *args, **kwargs ):
        super(PyTool, self).__init__( *args, **kwargs )
        if not hasattr( self, "func" ):
            raise Exception("A PyTool needs a 'func' attribute")
    def _run( self ):
        self.func()


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
        run_command( cmd, log=log_file, verbose=self.verbose )




class ScriptMixin( object ):
    def _make_script_file( self, **values_dict ):
        tmpl_file = os.path.join( self.tmpl_dir, self.tmpl_file )
        with open( tmpl_file, "r" ) as fp:
            tmpl_str = fp.read()
        script_file = os.path.join( self.output_dir, self.tmpl_file )
        with open( script_file, "w" ) as fp:
            fp.write( Template( tmpl_str ).substitute( **values_dict ) )
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






