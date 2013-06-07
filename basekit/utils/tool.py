#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division



import sys
import os
import argparse
import functools
import operator
import itertools
import inspect
import json
import cPickle as pickle
import csv
import string
import collections
import logging
import signal
import multiprocessing

import basekit.utils as utils
from basekit.utils import (
    try_int, get_index, boolean, working_directory, dir_walker, copy_dict,
    DefaultOrderedDict
)
from basekit.utils.timer import Timer
from basekit.utils.job import run_command
from basekit.utils.db import get_pdb_files


TIMEOUT_CMD = "timeout"


logging.basicConfig()
LOG = logging.getLogger('tool')
# LOG.setLevel( logging.ERROR )
LOG.setLevel( logging.WARNING )
# LOG.setLevel( logging.DEBUG )



class ToolParser( argparse.ArgumentParser ):
    def __init__( self, tool_class=None, description="", **kwargs ):
        if tool_class:
            description = tool_class.__doc__
            kwargs[ "add_help" ] = False
        else:
            if description:
                description += " The tools are accessible by the subcommands given below."
            else:
                description = "A collection of tools, accessible by the subcommands given below."
        super(ToolParser, self).__init__(
            description=description, epilog="basekit", **kwargs
        )
        if tool_class:
            group = self.add_argument_group( title="general arguments" )
            if not tool_class.no_output:
                group.add_argument( '-o', dest="output_dir", metavar='OUTPUT_DIR', type=str, default="./" )
                group.add_argument( '-t', dest="timeout", metavar='TIMEOUT', type=int, default=0 )
                group.add_argument( '-v', '--verbose', action='store_true' )
                group.add_argument( '-c', '--check', action='store_true' )
                group.add_argument( '-f', '--fileargs', action='store_true' )
            group.add_argument( '-h', '--help', action="help", help="show this help message and exit" )
    def error(self, message):
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def get_argument( params ):
    kwargs = {
        "default": params.get( "default_value", None ),
        "help": params.get( "help" )
    }
    if params.get('fixed'):
        kwargs["type"] = float
    elif params["type"]=="slider":
        kwargs["type"] = int
    elif params["type"] in [ "file", "text", "select" ]:
        kwargs["type"] = str
    elif params["type"]=="checkbox":
        if kwargs["default"]==False:
            kwargs["action"] = "store_true"
        else:
            kwargs["type"] = boolean
    return kwargs


def make_parser( Tool, parser=None ):
    if not parser:
        parser = ToolParser( description=Tool.__doc__ )
    arg_groups = DefaultOrderedDict( collections.OrderedDict )
    for name, params in Tool.args.iteritems():
        arg_groups[ params.get( "group" ) ][ name ] = params
    for group_name, args in arg_groups.iteritems():
        if group_name:
            group = parser.add_argument_group( title="%s arguments" % group_name )
        else:
            group = parser
        for name, params in args.iteritems():
            option = '--%s'%name if "default_value" in params else name
            group.add_argument( option, **get_argument( params ) )
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
    parser = ToolParser( description=description )
    subparsers = parser.add_subparsers( title='subcommands' )
    for name, Tool in tools.iteritems():
        subp = subparsers.add_parser( name, tool_class=Tool )
        make_parser( Tool, subp )
        subp.set_defaults(Tool=Tool)
    pargs = vars( parser.parse_args() )
    Tool = pargs.pop("Tool")
    args, kwargs = parse_args( Tool, kwargs=pargs )
    return Tool, args, kwargs






class Mixin( object ):
    pass


class TmplMixin( Mixin ):
    def _make_file_from_tmpl( self, tmpl_name, **values_dict ):
        tmpl_file = os.path.join( self.tmpl_dir, tmpl_name )
        with open( tmpl_file, "r" ) as fp:
            tmpl_str = fp.read()
        out_file = os.path.join( self.output_dir, tmpl_name )
        with open( out_file, "w" ) as fp:
            fp.write( string.Template( tmpl_str ).substitute( **values_dict ) )
        return out_file


class ScriptMixin( TmplMixin ):
    def _make_script_file( self, **values_dict ):
        return self._make_file_from_tmpl( self.script_tmpl, **values_dict )


class ProviMixin( TmplMixin ):
    def _make_provi_file( self, provi_tmpl=None, **values_dict ):
        provi_tmpl = provi_tmpl or self.provi_tmpl
        return self._make_file_from_tmpl( provi_tmpl, **values_dict )


class RecordsMixin( Mixin ):
    args = [
        { "name": "out_type", "type": "select", "options": [ "json", "csv", "pickle" ], 
          "default_value": "json" }
    ]
    def _init_records( self, stem, out_type="json", **kwargs ):
        if not hasattr( self, "RecordsClass" ):
            raise Exception("A RecordsMixin needs a 'RecordsClass' attribute")
        self.records = None
        self.out_type = out_type
        stem = stem or "%s_records" % self.name
        out_var = "%s_file" % self.out_type
        self.__dict__[ out_var ] = self.outpath( "%s.%s" % (stem, self.out_type) )
        self.output_files.append( self.__dict__[ out_var ] )
        if self.fileargs:
            self.read()
    def write_csv( self ):
        with open( self.csv_file, "w" ) as fp:
            cw = csv.writer( fp, delimiter=',')
            for r in self.records:
                cw.writerow( r )
    def write_json( self ):
        records_list = map( operator.methodcaller( "_asdict" ), self.records )
        with open( self.json_file, "w" ) as fp:
            json.dump( records_list, fp, indent=4 )
    def write_pickle( self ):
        with open( self.pickle_file, "w" ) as fp:
            pickle.dump( self.records, fp )
    def write( self ):
        getattr( self, "write_%s" % self.out_type )()
    def read_csv( self ):
        with open( self.csv_file, "r" ) as fp:
            self.records = map( self.RecordsClass._make, csv.reader( fp, delimiter=',') )
    def read_json( self ):
        with open( self.json_file, "r" ) as fp:
            records_list = json.load( fp, object_pairs_hook=collections.OrderedDict )
        self.records = map( lambda x: self.RecordsClass._make( x.itervalues() ), records_list )
    def read_pickle( self ):
        with open( self.pickle_file, "r" ) as fp:
            self.records = pickle.load( fp )
    def read( self ):
        getattr( self, "read_%s" % self.out_type )()


def call( tool ):
    try:
        return tool()
    except Exception as e:
        LOG.error( "[%s] %s" % ( tool.id, e ) )
    return tool


class ParallelMixin( Mixin ):
    args = [
        { "name": "parallel", "type": "select", "default_value": False,
          "options": [ "directory", "pdb_archive", "list" ] },
        { "name": "sample", "type": "slider", "range": [0, 100], "default_value": None },
        { "name": "sample_start", "type": "slider", "range": [0, 100], "default_value": 0 }
    ]
    def _init_parallel( self, file_path, parallel=None, sample=None, sample_start=0, **kwargs ):
        self.file_path = self.abspath( file_path )
        self.parallel = parallel
        self.sample_start = sample_start or 0
        self.sample_end = None if not sample else self.sample_start+sample
    def _make_tool_list( self ):
        if self.parallel=="pdb_archive":
            file_list = get_pdb_files( self.file_path, pattern=".pdb" )
        elif self.parallel in [ "directory", "dir" ]:
            file_list = map( operator.itemgetter(1), dir_walker( self.file_path, pattern=".+\.pdb" ) )
        elif self.parallel=="list":
            file_list = self.file_path.split()
        else:
            raise Exception( "unknown value '%s' for 'parallel'" % self.parallel )
        file_list = itertools.islice( file_list, self.sample_start, self.sample_end )
        tool_list = []
        kwargs = { 
            "run": False
        }
        for input_file in file_list:
            stem = utils.path.stem( input_file )
            output_dir = self.outpath( os.path.join( "parallel", stem ) )
            tool_list.append( self.ParallelClass(
                input_file, pdb_id=stem, **copy_dict( kwargs, output_dir=output_dir )
            ))
        self.tool_list = tool_list
    def _func_parallel( self, nworkers=None ):
        # !important - allows one to abort via CTRL-C
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        multiprocessing.log_to_stderr( logging.ERROR )
        
        if not nworkers: nworkers = multiprocessing.cpu_count()
        p = multiprocessing.Pool( nworkers, maxtasksperchild=50 )

        data = p.imap( call, self.tool_list )

        p.close()
        p.join()

        return data




class ToolMetaclass(type):
    def __init__(cls, name, bases, dct):
        # print cls, name, bases, dct
        if not "no_output" in dct:
            cls.no_output = False

        args = collections.OrderedDict()
        for a in dct.get( "args", [] ):
            args[ a.pop("name") ] = a

        if RecordsMixin in bases:
            for a in RecordsMixin.__dict__.get( "args", [] ):
                a = a.copy()
                a["group"] = "records"
                args[ a.pop("name") ] = a

        if ParallelMixin in bases:
            for a in ParallelMixin.__dict__.get( "args", [] ):
                a = a.copy()
                a["group"] = "parallel"
                args[ a.pop("name") ] = a
            if not "ParallelClass" in dct:
                cls.ParallelClass = cls

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
        
        if kwargs.get("run", True) and not kwargs.get("check", False) and not self.fileargs:
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


class DbTool( Tool ):
    def __init__( self, cmd=None, *args, **kwargs ):
        super(DbTool, self).__init__( *args, **kwargs )
        if not hasattr( self, "schema" ):
            raise Exception("A DbTool needs a 'schema' attribute")
    def _run( self ):
        pass










