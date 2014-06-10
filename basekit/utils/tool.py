from __future__ import with_statement
from __future__ import division



import sys
import os
import re
import argparse
import operator
import itertools
import cStringIO
import json
import cPickle as pickle
import csv
import sqlite3
import string
import collections
import logging
import signal
import multiprocessing

import jinja2

import basekit.utils.numpdb as numpdb
from basekit import utils
from basekit.utils import (
    boolean, working_directory, dir_walker, copy_dict,
    DefaultOrderedDict
)
from basekit.utils.job import run_command
from basekit.utils.db import get_pdb_files


def _dir_init( tool_path, tool_name ):
    """use like this:
        DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "dowser" )
    """
    DIR = os.path.split( os.path.abspath(tool_path) )[0]
    PARENT_DIR = os.path.split( DIR )[0]
    TMPL_DIR = os.path.join( PARENT_DIR, "data", tool_name )
    return DIR, PARENT_DIR, TMPL_DIR


TIMEOUT_CMD = "timeout"


logging.basicConfig()
LOG = logging.getLogger('tool')
LOG.setLevel( logging.INFO )



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
                    type=str, default="" 
                )
                group.add_argument( 
                    '-t', dest="timeout", metavar='TIMEOUT', 
                    type=int, default=0 
                )
                group.add_argument( '-v', '--verbose', action='store_true' )
                group.add_argument( '-c', '--check', action='store_true' )
                group.add_argument( '-a', '--fileargs', action='store_true' )
                group.add_argument( '-d', '--debug', action='store_true' )
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
    if "options" in params:
        kwargs["choices"] = params["options"]

    if params.get( "nargs" ):
        kwargs["nargs"] = params["nargs"]

    if params["type"] in [ "float" ]:
        kwargs["type"] = float
    elif params["type"] in [ "int", "integer" ]:
        kwargs["type"] = int
    elif params["type"] in [ "file", "dir", "sele", "str", "string", "list" ]:
        kwargs["type"] = str
    elif params["type"] in [ "bool", "boolean" ]:
        if kwargs["default"]==False:
            kwargs["action"] = "store_true"
        else:
            kwargs["type"] = boolean
    else:
        raise Exception( "type '%s' not known" % params["type"] )
    
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

        def make_out( params ):
            return params
        out = collections.OrderedDict()
        for p in dct.get( "out", [] ):
            out[ p["name"] ] = make_out( p )

        for mixin_name, mixin_cls in MIXIN_REGISTER.iteritems():
            if mixin_cls in bases:
                for p in mixin_cls.__dict__.get( "args", [] ):
                    p["group"] = mixin_cls.__name__[:-5].lower()
                    p = make_arg( p )
                    args[ p["name"] ] = p
                for p in mixin_cls.__dict__.get( "out", [] ):
                    out[ p["name"] ] = make_out( p )

        cls.args = args
        cls.out = out

        if MIXIN_REGISTER['ParallelMixin'] in bases:
            if not "ParallelClass" in dct:
                cls.ParallelClass = cls
            cls_run = cls._run
            cls._run = lambda self: MIXIN_REGISTER['ParallelMixin']._run(
                self, cls_run
            )

        if not "no_output" in dct:
            cls.no_output = False


class Mixin( object ):
    __metaclass__ = ToolMetaclass


class TmplMixin( Mixin ):
    def _make_file_from_tmpl( self, tmpl_name, output_dir=None, out_name=None,
                              prefix=None, use_jinja2=False, **values_dict ):
        output_dir = output_dir or self.output_dir
        tmpl_file = os.path.join( self.tmpl_dir, tmpl_name )
        with open( tmpl_file, "r" ) as fp:
            tmpl_str = fp.read()
        out_file = os.path.join( output_dir, tmpl_name )
        out_file = utils.path.mod( out_file, prefix=prefix, name=out_name )
        if use_jinja2:
            out = jinja2.Template( tmpl_str ).render( **values_dict )
        else:
            out = string.Template( tmpl_str ).substitute( **values_dict )
        with open( out_file, "w" ) as fp:
            fp.write( out )
        return out_file


class ScriptMixin( TmplMixin ):
    def _make_script_file( self, **kwargs ):
        return self._make_file_from_tmpl( self.script_tmpl, **kwargs )


class ProviMixin( TmplMixin ):
    def _make_provi_file( self, provi_tmpl=None, **kwargs ):
        provi_tmpl = provi_tmpl or self.provi_tmpl
        return self._make_file_from_tmpl( provi_tmpl, **kwargs )



BACKEND_REGISTER = {}
class BackendMetaclass( type ):
    def __init__(cls, name, bases, dct):
        BACKEND_REGISTER[ cls.name ] = cls

class RecordsBackend( object ):
    __metaclass__ = BackendMetaclass
    name = "dat"
    def __init__( self, file_name, cls ):
        self.file_name = file_name
        self.cls = cls
    def write( self, records ):
        pass
    def read( self ):
        pass 

class CsvBackend( RecordsBackend ):
    name = "csv"
    def write( self, records ):
        with open( self.file_name, "w" ) as fp:
            cw = csv.writer( fp, delimiter=',')
            cw.writerow( self.cls._fields )
            for r in records:
                cw.writerow( r )
    def read( self ):
        with open( self.file_name, "r" ) as fp:
            cr = csv.reader( fp, delimiter=',')
            cr.next() # ignore header
            return map( self.cls._make, cr )

class JsonBackend( RecordsBackend ):
    name = "json"
    def write( self, records ):
        records_list = map( 
            operator.methodcaller( "_asdict" ), records
        )
        with open( self.file_name, "w" ) as fp:
            json.dump( records_list, fp, indent=4 )
    def read( self ):
        with open( self.file_name, "r" ) as fp:
            records_list = json.load( 
                fp, object_pairs_hook=collections.OrderedDict
            )
        return map( 
            lambda x: self.cls._make( x.itervalues() ), 
            records_list 
        ) 

class PickleBackend( RecordsBackend ):
    name = "pickle"
    def write( self, records ):
        with open( self.file_name, "w" ) as fp:
            pickle.dump( records, fp )
    def read( self ):
        with open( self.file_name, "r" ) as fp:
            return pickle.load( fp ) 

class SqliteBackend( RecordsBackend ):
    name = "sqlite"
    conn_kwargs = { "isolation_level": "EXCLUSIVE" }
    def write( self, records ):
        utils.path.remove( self.file_name )
        qn = ",".join( "?" * len( self.cls._fields ) )
        with sqlite3.connect( self.file_name, **self.conn_kwargs ) as conn:
            c = conn.cursor()
            c.execute( 
                'CREATE TABLE %s (%s)' % ( 
                    self.cls.__name__, ",".join( self.cls._fields ) 
                )
            )
            c.executemany( 
                'INSERT INTO %s VALUES (%s)' % ( self.cls.__name__, qn ),
                itertools.imap( tuple, records ) 
            )
    def read( self ):
        return self.query()
    def query( self, where=None, count=False, limit=None ):
        self.q = 'SELECT\n\t%s\nFROM\n\t%s' % (
            'COUNT(*)' if count else '*', self.cls.__name__
        )
        if where:
            self.q += '\nWHERE\n\t%s' % where
        if limit:
            self.q += '\nLIMIT\n\t%i' % limit
        with sqlite3.connect( self.file_name, **self.conn_kwargs ) as conn:
            if not count:
                conn.row_factory = lambda x, y: self.cls( *y )
            c = conn.cursor()
            c.execute( self.q )
            return c.fetchone() if count else c.fetchall()

def records_backend( backend, file_name, records_cls ):
    backend_cls = BACKEND_REGISTER[ backend ]
    return backend_cls( file_name, records_cls )

class RecordsMixin( Mixin ):
    args = [
        _( "backend|rb", type="str", default="json", 
            options=BACKEND_REGISTER.keys(), metavar="B" ),
        # TODO only when ParallelMixin
        _( "backend_parallel|rbp", type="str", default="json", 
            options=BACKEND_REGISTER.keys(), metavar="B" )
    ]
    def _init_records( self, input_file, **kwargs ):
        if not hasattr( self, "RecordsClass" ):
            raise Exception("A RecordsMixin needs a 'RecordsClass' attribute")
        self.records = []
        if input_file:
            stem = utils.path.stem( input_file ) 
        else:
            stem = "%s_records" % self.name
        self.backend_obj = self._backend( stem, self.backend )
        # TODO only when ParallelMixin
        self.backend_obj_p = self._backend( stem, self.backend_parallel )
        if self.fileargs:
            self.read()
    def _backend( self, stem, backend ):
        backend_cls = BACKEND_REGISTER[ backend ]
        file_name = os.path.join( 
            self.output_dir, "%s.%s" % ( stem, backend_cls.name )
        )
        return backend_cls( file_name, self.RecordsClass )
    def write( self ):
        self.backend_obj.write( self.records )
    def read( self ):
        try:
            self.records = self.backend_obj.read()
        except Exception as e:
            print e
    # TODO only when ParallelMixin
    def _parallel_results( self, tool_list ):
        self.records = list(itertools.chain.from_iterable(
            map( operator.attrgetter( "records" ), tool_list )
        ))
        self.backend_obj_p.write( self.records )



def call( tool ):
    try:
        tool()
        if hasattr( tool, "check_only" ) and tool.check_only:
            info = tool.check( full=True )
        else:
            info = str( tool )
        LOG.info( "[%s] %s" % ( tool.id, info ) )
    except Exception as e:
        LOG.error( "[%s] %s" % ( tool.id, e ) )
        if tool.debug:
            import traceback
            traceback.print_exc()
    return tool


class ParallelMixin( Mixin ):
    args = [
        _( "parallel|p", type="str", default=False,
            options=[ "", "directory", "pdb_archive", "list", "file", "file_list", "data" ] ),
        _( "interval|i", type="int", default=[ 0, None ], 
            metavar=("BEG", "END"), nargs=2 ),
        _( "filter_id|id", type="str", default=[], nargs="*" ),
        _( "nworkers|nw", type="int", default=0 ),
    ]
    def _init_parallel( self, input_data, **kwargs ):
        if not hasattr( self, "_parallel_results" ):
            raise Exception( "'_parallel_results' method required" )
        if self.parallel in [ "pdb_archive", "directory", "dir", "file" ]:
            self.input_data = self.abspath( input_data )
        elif self.parallel in [ "file_list" ]:
            self.input_data = map(
                self.abspath,
                re.split( "[\s,]+", input_data.strip() )
            )
        elif self.parallel in [ "list" ]:
            self.input_data = re.split( "[\s,]+", input_data.strip() )
        else:
            self.input_data = input_data
        self.tool_kwargs = kwargs
        if self.fileargs and self.parallel:
            self.tool_kwargs[ "fileargs" ] = True
            self._make_tool_list()
            self._parallel_results( self.tool_list )
        if not self.nworkers:
            self.nworkers = multiprocessing.cpu_count()
    def _make_tool_list( self ):
        if self.parallel=="pdb_archive":
            input_list = get_pdb_files( 
                self.input_data, pattern=".pdb" )
        elif self.parallel in [ "directory", "dir" ]:
            input_list = map( 
                operator.itemgetter(1), 
                dir_walker( self.input_data, pattern=".+\.pdb" ) 
            )
        elif self.parallel in [ "file" ]:
            input_list = []
            with open( self.input_data, "r" ) as fp:
                for line in fp:
                    ls = line.strip()
                    if ls: 
                        input_list.append( ls )
        elif self.parallel in [ "list", "file_list", "data" ]:
            input_list = self.input_data
        else:
            raise Exception( 
                "unknown value '%s' for 'parallel'" % self.parallel 
            )
        input_list = itertools.islice( 
            input_list, self.interval[0], self.interval[1]
        )
        tool_list = []
        self.tool_kwargs["run"] = False
        for i, input_elm in enumerate(input_list):
            try:
                stem = utils.path.stem( input_elm )
            except:
                stem = str(i)
            output_dir = self.outpath( os.path.join( "parallel", stem ) )
            tool = self.ParallelClass( input_elm, **copy_dict( 
                self.tool_kwargs, parallel=False,
                output_dir=output_dir,
            ))
            tool.id = stem
            tool_list.append( tool )
        if self.filter_id:
            tool_list = filter( lambda x: x.id in self.filter_id, tool_list )
        self.tool_list = tool_list
    def _func_parallel( self ):
        # !important - allows one to abort via CTRL-C
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        multiprocessing.log_to_stderr( logging.ERROR )
        p = multiprocessing.Pool( self.nworkers, maxtasksperchild=50 )
        data = p.imap( call, self.tool_list )
        p.close()
        p.join()
        return list( data )
    def _run( self, fn ):
        if self.parallel:
            self._make_tool_list()
            self.tool_results = self._func_parallel()
            self._parallel_results( self.tool_results )
        else:
            fn( self )



class Tool( object ):
    __metaclass__ = ToolMetaclass
    def __init__( self, *args, **kwargs ):
        args, kwargs = self._pre_init( args, kwargs )
        self.name = self.__class__.__name__.lower()

        # hidden kwargs
        self.pre_exec = kwargs.get("pre_exec", True)
        self.post_exec = kwargs.get("post_exec", True)

        # general
        self.timeout = kwargs.get("timeout", None)
        self.fileargs = kwargs.get("fileargs", False)
        self.verbose = kwargs.get("verbose", False)
        self.debug = kwargs.get("debug", False)
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

        self.params = { "args": args, "kwargs": kwargs }

        self.input_files_dict = { "cls": self }
        args_iter = iter(args)
        for name, params in self.args.iteritems():
            if "default" in params:
                value = kwargs.get( name, params["default"] )
            else:
                try:
                    value = args_iter.next()
                except StopIteration:
                    raise Exception( 'Missing positional argument' )
            # TODO check if the name already exists
            self.__dict__[ name ] = self.__prep_arg( value, params )
            if params["type"]=="file" and "nargs" not in params:
                self.input_files_dict[ name ] = utils.Bunch(
                    stem=utils.path.stem( value ),
                    basename=os.path.basename( value ),
                    ext=utils.path.ext( value )
                )
            elif params["type"]=="str" and "nargs" not in params:
                self.input_files_dict[ name ] = value
        
        self.output_files = []
        for name, params in self.out.iteritems():
            value = self.__prep_out( params )
            # TODO check if the name already exists
            self.__dict__[ name ] = value
            if "file" in params and not params.get( "optional", False ):
                self.output_files.append( value )

        if hasattr( self, "stdout_file" ):
            raise Exception(
                "'stdout_file' attribute already defined (%s)" % self.name
            )
        self.stdout_file = self.outpath( "%s_cmd.log" % self.name )

        self._init( *args, **kwargs )
        
        if kwargs.get("run", True) and not kwargs.get("check", False) and not self.fileargs:
            self.__run()
    def __prep_arg( self, value, params ):
        if not value: 
            return value
        if params.get("type") in [ "file", "dir" ]:
            if "nargs" in params:
                return map( self.abspath, value )
            return self.abspath( value )
        elif params.get("type")=="sele":
            return numpdb.numsele( value )
        return value
    def __prep_out( self, params ):
        if "file" in params:
            return self.outpath( params["file"] )
        elif "dir" in params:
            return self.subdir( params["dir"] )
        raise Exception( "do not know how to prep output" )
    def __run( self ):
        with working_directory( self.output_dir ):
            if self.pre_exec: self._pre_exec()
            self._run()
            if self.post_exec: self._post_exec()
    def __call__( self ):
        self.__run()
        return self
    def _pre_init( self, args, kwargs ):
        return args, kwargs
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
            return os.path.isfile(f)# and os.path.getsize(f)>0
        except:
            return False
    def check( self, full=False ):
        if hasattr( self, "output_files" ):
            with working_directory( self.output_dir ):
                if full:
                    for of in self.output_files:
                        if not self.__check_file( of ):
                            return self.relpath( of )
                    else:
                        return "Ok"
                else:
                    isfile = map( self.__check_file, self.output_files )
                    return all( isfile )
        else:
            return True
    def __str__( self ):
        status = "ok" if self.check() else "failed"
        return "%s status: %s" % ( self.name, status )
    def relpath( self, path, no_ext=False ):
        if no_ext:
            path = os.path.splitext( path )[0]
        return str( os.path.relpath( path, self.output_dir ) )
    def outpath( self, file_name ):
        file_name = file_name.format( **self.input_files_dict )
        return os.path.join( self.output_dir, file_name )
    def abspath( self, path ):
        return os.path.abspath( path )
    def subdir( self, directory, filename=None ):
        subdir = os.path.join( self.output_dir, directory )
        if not os.path.exists( subdir ):
            os.makedirs( subdir )
        if filename:
            return os.path.join( subdir, filename )
        else:
            return subdir
    def datapath( self, file_name ):
        if not hasattr( self, "tmpl_dir" ):
            raise "No data path available."
        return os.path.join( self.tmpl_dir, file_name )




class PyTool( Tool ):
    def __init__( self, *args, **kwargs ):
        super(PyTool, self).__init__( *args, **kwargs )
        if not hasattr( self, "func" ):
            raise Exception("A PyTool needs a 'func' attribute")
    def _run( self ):
        if self.verbose:
            return self.func()
        else:
            stdout, stderr = sys.stdout, sys.stderr
            log = cStringIO.StringIO()
            sys.stdout = log
            sys.stderr = log
            try:
                return self.func()
            finally:
                sys.stdout, sys.stderr = stdout, stderr


class CmdTool( Tool ):
    cmd_input = None
    use_shell = False
    no_cmd = None
    def __init__( self, *args, **kwargs ):
        super(CmdTool, self).__init__( *args, **kwargs )
        if not hasattr( self, "cmd" ):
            raise Exception("A CmdTool needs a 'cmd' attribute")
    def _run( self ):
        if self.no_cmd:
            return
        cmd = self.cmd
        if self.verbose:
            print " ".join( map( str, cmd ) )
        if self.timeout:
            cmd = [ TIMEOUT_CMD, self.timeout ] + cmd
        ret = run_command( 
            cmd, log=self.stdout_file, verbose=self.verbose, 
            input_file=self.cmd_input, shell=self.use_shell
        )
        return ret










