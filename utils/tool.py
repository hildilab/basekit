#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division



import re
import os
import shutil
import argparse
import sqlite3
import functools
import itertools
from string import Template

from utils import try_int, get_index, boolean
from utils.timer import Timer
from utils.job import run_command2, working_directory


TIMEOUT_CMD = "timeout"


class Tool( object ):
    def __init__( self, output_dir=".", run=True, timeout=None ):
        output_dir = os.path.abspath( output_dir )
        if not os.path.exists( output_dir ):
            os.makedirs( output_dir )
        self.output_dir = output_dir
        self.timeout = timeout
        if run:
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
    def __init__( self, **kw ):
        if not hasattr( self, "func" ):
            raise Exception("A PyTool needs a 'func' attribute")
        super(PyTool, self).__init__( **kw )
    def _run( self ):
        self.func()


class CmdTool( Tool ):
    def __init__( self, **kw ):
        if not hasattr( self, "cmd" ):
            raise Exception("A CmdTool needs a 'cmd' attribute")
        super(CmdTool, self).__init__( **kw )
    def _run( self ):
        cmd = self.cmd
        if self.timeout:
            cmd = [ TIMEOUT_CMD, self.timeout ] + cmd
        log_file = "%s_cmd.log" % self.__class__.__name__.lower()
        run_command2( cmd, log=log_file )




