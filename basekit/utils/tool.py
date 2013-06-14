#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division



import sys
import os
import re
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

import basekit.utils.numpdb as numpdb
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




def _( name, **kwargs ):
    kwargs.update( name=name )
    return kwargs


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
                group.add_argument( 
                    '-o', dest="output_dir", metavar='OUTPUT_DIR', 
                    type=str, default="./" 
                )
                group.add_argument( 
                    '-t', dest="timeout", metavar='TIMEOUT', 
                    type=int, default=0 
                )
                group.add_argument( '-v', '--verbose', action='store_true' )
                group.add_argument( '-c', '--check', action='store_true' )
                group.add_argument( '-a', '--fileargs', action='store_true' )
            group.add_argument( 
                '-h', '--help', action="help", 
                help="show this help message and exit" 
            )
    def error(self, message):
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def get_argument( params ):
    kwargs = {
        "default": params.get( "default" ),
        "help": params.get( "help" )
    }

    if "required" in params:
        kwargs["required"] = params["required"]
    if "metavar" in params:
        kwargs["metavar"] = params["metavar"]
    if "dest" in params:
        kwargs["dest"] = params["dest"]
    if "action" in params:
        kwargs["action"] = params["action"]
    if "choices" in params:
        kwargs["choices"] = params["choices"]

    if params.get( "nargs" ):
        kwargs["nargs"] = params["nargs"]

    if params.get('fixed') or params["type"] in [ "float" ]:
        kwargs["type"] = float
    elif params["type"] in [ "slider", "int" ]:
        kwargs["type"] = int
    elif params["type"] in [ "file", "dir", "text", "select", "str" ]:
        kwargs["type"] = str
    elif params["type"]=="checkbox":
        if kwargs["default"]==False:
            kwargs["action"] = "store_true"
        else:
            kwargs["type"] = boolean
    
    return kwargs


def make_parser( Tool, parser=None ):
    if not parser:
        parser = ToolParser( tool_class=Tool, description=Tool.__doc__ )
    arg_groups = DefaultOrderedDict( collections.OrderedDict )
    for name, params in Tool.args.iteritems():
        arg_groups[ params.get( "group" ) ][ name ] = params
    for group_name, args in arg_groups.iteritems():
        if group_name:
            group = parser.add_argument_group( 
                title="%s arguments" % group_name 
            )
        else:
            group = parser
        for name, params in args.iteritems():
            group.add_argument( 
                *params["flags"], **get_argument( params ) 
            )
    return parser


def parse_args( Tool, kwargs=None ):
    if not kwargs:
        parser = make_parser( Tool )
        kwargs = vars( parser.parse_args() )
    args = []
    for name, params in Tool.args.iteritems():
        if "default" not in params:
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


MIXIN_REGISTER = {}
class ToolMetaclass( type ):
    def __init__(cls, name, bases, dct):
        # print cls, name, bases, dct
        if name.endswith("Mixin"):
            MIXIN_REGISTER[ name ] = cls
            return

        def make_arg( params ):
            p = params.copy()
            flags = p["name"].split("|")
            p["name"] = flags[0]
            if "default" in p:
                flags[0] = "--%s" % flags[0]
            if len(flags)>1:
                if "default" not in p:
                    flags[0] = "--%s" % flags[0]
                    p["required"] = True
                for i in range( 1, len(flags) ):
                    flags[i] = "-%s" % flags[i]
            p["flags"] = flags
            return p

        args = collections.OrderedDict()
        for p in dct.get( "args", [] ):
            p = make_arg( p )
            args[ p["name"] ] = p

        for mixin_name, mixin_cls in MIXIN_REGISTER.iteritems():
            if mixin_cls in bases:
                for p in mixin_cls.__dict__.get( "args", [] ):
                    p["group"] = mixin_cls.__name__[:-5].lower()
                    p = make_arg( p )
                    args[ p["name"] ] = p

        cls.args = args

        if MIXIN_REGISTER['ParallelMixin'] in bases:
            if not "ParallelClass" in dct:
                cls.ParallelClass = cls

        # TODO remove
        if not "no_output" in dct:
            cls.no_output = False

        def make_out( params ):
            return params

        out = collections.OrderedDict()
        for p in dct.get( "out", [] ):
            out[ p["name"] ] = make_out( p )

        cls.out = out



class Mixin( object ):
    __metaclass__ = ToolMetaclass


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
        { "name": "records_type", "type": "select", 
            "options": [ "json", "csv", "pickle" ], 
            "default": "json" }
    ]
    def _init_records( self, input_file, **kwargs ):
        if not hasattr( self, "RecordsClass" ):
            raise Exception("A RecordsMixin needs a 'RecordsClass' attribute")
        self.records = None
        if input_file:
            stem = utils.path.stem( input_file ) 
        else:
            stem = "%s_records" % self.name
        out_var = "records_%s" % self.records_type
        self.__dict__[ out_var ] = self.outpath( 
            "%s.%s" % ( stem, self.records_type ) 
        )
        self.output_files.append( self.__dict__[ out_var ] )
        if self.fileargs:
            self.read()
    def write_csv( self ):
        with open( self.records_csv, "w" ) as fp:
            cw = csv.writer( fp, delimiter=',')
            for r in self.records:
                cw.writerow( r )
    def write_json( self ):
        records_list = map( 
            operator.methodcaller( "_asdict" ), 
            self.records 
        )
        with open( self.records_json, "w" ) as fp:
            json.dump( records_list, fp, indent=4 )
    def write_pickle( self ):
        with open( self.records_pickle, "w" ) as fp:
            pickle.dump( self.records, fp )
    def write( self ):
        getattr( self, "write_%s" % self.records_type )()
    def read_csv( self ):
        with open( self.records_csv, "r" ) as fp:
            self.records = map( 
                self.RecordsClass._make, 
                csv.reader( fp, delimiter=',') 
            )
    def read_json( self ):
        with open( self.records_json, "r" ) as fp:
            records_list = json.load( 
                fp, object_pairs_hook=collections.OrderedDict
            )
        self.records = map( 
            lambda x: self.RecordsClass._make( x.itervalues() ), 
            records_list 
        )
    def read_pickle( self ):
        with open( self.records_pickle, "r" ) as fp:
            self.records = pickle.load( fp )
    def read( self ):
        getattr( self, "read_%s" % self.records_type )()



def call( tool ):
    try:
        return tool()
    except Exception as e:
        LOG.error( "[%s] %s" % ( tool.id, e ) )
    return tool


class ParallelMixin( Mixin ):
    args = [
        { "name": "parallel", "type": "select", "default": False,
          "options": [ "directory", "pdb_archive", "list" ] },
        { "name": "sample", "type": "slider", "range": [0, 100], 
            "default": None },
        { "name": "sample_start", "type": "slider", "range": [0, 100], 
            "default": 0 }
    ]
    def _init_parallel( self, file_input, parallel=None, 
                        sample=None, sample_start=0, **kwargs ):
        if parallel in [ "pdb_archive", "directory", "dir" ]:
            self.file_input = self.abspath( file_input )
        elif parallel in [ "list" ]:
            self.file_input = map(
                self.abspath,
                re.split( "[\s,]+", file_input.strip() )
            )
        else:
            self.file_input = file_input
        self.parallel = parallel
        self.sample_start = sample_start or 0
        self.sample_end = None if not sample else self.sample_start+sample
        self.tool_kwargs = kwargs
    def _make_tool_list( self ):
        if self.parallel=="pdb_archive":
            file_list = get_pdb_files( 
                self.file_input, pattern=".pdb" )
        elif self.parallel in [ "directory", "dir" ]:
            file_list = map( 
                operator.itemgetter(1), 
                dir_walker( self.file_input, pattern=".+\.pdb" ) 
            )
        elif self.parallel=="list":
            file_list = self.file_input
        else:
            raise Exception( 
                "unknown value '%s' for 'parallel'" % self.parallel 
            )
        file_list = itertools.islice( 
            file_list, self.sample_start, self.sample_end
        )
        tool_list = []
        self.tool_kwargs["run"] = False
        for input_file in file_list:
            stem = utils.path.stem( input_file )
            output_dir = self.outpath( os.path.join( "parallel", stem ) )
            tool_list.append( self.ParallelClass(
                input_file, **copy_dict( 
                    self.tool_kwargs, output_dir=output_dir,
                )
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






class Tool( object ):
    __metaclass__ = ToolMetaclass
    def __init__( self, *args, **kwargs ):
        self.name = self.__class__.__name__.lower()

        self.input_files_dict = { "cls": self }
        args_iter = iter(args)
        for name, params in self.args.iteritems():
            if "default" in params:
                value = kwargs.get( name, params["default"] )
            else:
                value = args_iter.next()
            # TODO check if the name already exists
            self.__dict__[ name ] = self.__prep_arg( value, params )
            if params["type"]=="file":
                self.input_files_dict[ name ] = utils.Bunch(
                    stem=utils.path.stem( value )
                )

        self.timeout = kwargs.get("timeout", None)
        self.fileargs = kwargs.get("fileargs", False)
        self.verbose = kwargs.get("verbose", False)
        self.output_dir = os.path.abspath( 
            kwargs.get("output_dir", ".") 
        ) + os.sep
        
        if not self.no_output:
            if not os.path.exists( self.output_dir ):
                os.makedirs( self.output_dir )
            self.args_file = os.path.join( 
                self.output_dir, "%s.json" % self.name 
            )
            if self.fileargs:
                with open( self.args_file, "r" ) as fp:
                    args, kwargs = json.load( fp )
            else:
                with open( self.args_file, "w" ) as fp:
                    json.dump( ( args, kwargs ), fp, indent=4 )
        
        self.output_files = []
        for name, params in self.out.iteritems():
            value = self.__prep_out( params )
            # TODO check if the name already exists
            self.__dict__[ name ] = value
            self.output_files.append( value )

        self._init( *args, **kwargs )
        
        if kwargs.get("run", True) and not kwargs.get("check", False) and not self.fileargs:
            self.__run()
    def __prep_arg( self, value, params ):
        if not value: 
            return value
        if params.get("type") in [ "file", "dir" ]:
            return self.abspath( value )
        elif params.get("type")=="sele":
            return numpdb.numsele( value )
        return value
    def __prep_out( self, params ):
        if "file" in params:
            return self.outpath( params["file"] )
        elif "dir" in params:
            return self.subdir( params["dir"] )
        raise "do not know how to prep output"
    def __run( self ):
        with working_directory( self.output_dir ):
            self._pre_exec()
            self._run()
            self._post_exec()
    def __call__( self ):
        self.__run()
        return self
    def _init( self, *args, **kwargs ):
        pass
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
        file_name = file_name.format( **self.input_files_dict )
        return os.path.join( self.output_dir, file_name )
    def abspath( self, path ):
        return os.path.abspath( path )
    def subdir( self, directory ):
        subdir = os.path.join( self.output_dir, directory )
        if not os.path.exists( subdir ):
            os.makedirs( subdir )
        return subdir




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










