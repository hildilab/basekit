from __future__ import division
from __future__ import with_statement

import multiprocessing
import functools
import os
import sys
import contextlib
import subprocess
import fcntl
import select
import time
import signal
import logging




@contextlib.contextmanager 
def working_directory(directory): 
    original_directory = os.getcwd()
    try: 
        os.chdir(directory) 
        yield directory 
    finally: 
        os.chdir(original_directory) 


# TODO: remove shell dependency
def run_command( program, log=None, stdout=None, close_fds=True ):
    # this function allow a stream-like access to stdout/stderr 
    command="(%s) 2>&1" % program
    p = subprocess.Popen( command, shell=True, stdout=subprocess.PIPE, close_fds=close_fds )
    if log: fp = open( log, 'w' )
    fcntl.fcntl(
        p.stdout.fileno(),
        fcntl.F_SETFL,
        fcntl.fcntl(p.stdout.fileno(), fcntl.F_GETFL) | os.O_NONBLOCK,
    )
    while True:
        readx = select.select([p.stdout.fileno()], [], [])[0]
        if readx:
            chunk = p.stdout.read()
            if chunk == '':
                break
            if stdout: sys.stdout.write( chunk )
            if log: fp.write( chunk )
        time.sleep(.1)
    if log: fp.close()
    return


def run_command2( cmd, cwd=".", log=None ):
    kwargs = {
        "args": cmd,
        "cwd": cwd,
        "stdout": os.devnull,
        "stderr": subprocess.STDOUT,
        "env": os.environ,
        "preexec_fn": os.setpgrp
    }
    if log:
        with open( log, "w" ) as fp:
            kwargs.update(stdout=fp)
            ret = subprocess.call( **kwargs )
    else:
        ret = subprocess.call( **kwargs )
    return ret


def _prep_func(fn):
    def wrapped( fpath ):
        try:
            return ( fpath, fn( fpath ) )
        except Exception, e:
            print e
            return ( fpath, None )
    return functools.update_wrapper( wrapped, fn )


def do_parallel( func, files, nworkers=None ):

    # !important - allows one to abort via CTRL-C
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    multiprocessing.log_to_stderr( logging.ERROR )
    
    if not nworkers:
        nworkers = multiprocessing.cpu_count()
    
    p = multiprocessing.Pool( nworkers )
    data = p.map( func, files )
    p.close()
    p.join()

    return data
