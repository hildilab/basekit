#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division



import sys
import re
import os
import shutil
import argparse
import sqlite3
import functools
import itertools
import inspect
from string import Template
from collections import OrderedDict

from basekit.utils import try_int, get_index, boolean
from basekit.utils.timer import Timer
from basekit.utils.job import run_command2, working_directory


TIMEOUT_CMD = "timeout"


def make_args( args ):
    _args = OrderedDict()
    for a in args:
        _args[ a.pop("name") ] = a
    return _args


class Tool( object ):
    _args = make_args([])
    def __init__( self, *args, **kwargs ):
        self.name = self.__class__.__name__.lower()

        self.timeout = kwargs.get("timeout", None)
        self.fileargs = kwargs.get("fileargs", False)

        self.output_dir = os.path.abspath( kwargs.get("output_dir", ".") )
        if not os.path.exists( self.output_dir ):
            os.makedirs( self.output_dir )

        self.args_file = os.path.join( self.output_dir, "%s.args" % self.name )
        if self.fileargs:
            with open( self.args_file, "r" ) as fp:
                args = fp.read().split(",")
        else:
            with open( self.args_file, "w" ) as fp:
                fp.write( ",".join(args) )

        self._init( *args )

        if kwargs.get("run", True):
            self.__run()
    def __run( self ):
        with working_directory( self.output_dir ):
            self._pre_exec()
            self._run()
    def __call__( self ):
        self.__run()
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
        run_command2( cmd, log=log_file )




