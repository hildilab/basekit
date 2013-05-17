from __future__ import with_statement

import multiprocessing
import functools
import os
import sys
import subprocess
import fcntl
import select
import time
import signal
import logging




def run_command( cmd, cwd=".", log=None, verbose=False ):
    kwargs = {
        "args": map( str, cmd ),
        "cwd": cwd,
        "stdout": os.devnull,
        "stderr": subprocess.STDOUT,
        "env": os.environ,
        "preexec_fn": os.setpgrp
    }
    if verbose:
        kwargs.update( stdout=subprocess.PIPE )
        if log:
            with open( log, "w" ) as fp:
                p = subprocess.Popen( **kwargs )
                fcntl.fcntl(
                    p.stdout.fileno(),
                    fcntl.F_SETFL,
                    fcntl.fcntl(p.stdout.fileno(), fcntl.F_GETFL) | os.O_NONBLOCK,
                )
                while True:
                    readx = select.select([p.stdout.fileno()], [], [])[0]
                    if readx:
                        chunk = p.stdout.read()
                        if chunk == '': break
                        sys.stdout.write( chunk )
                        fp.write( chunk )
                    time.sleep(.1)
            return p.returncode
        else:
            return subprocess.call( **kwargs )
    elif log:
        with open( log, "w" ) as fp:
            kwargs.update( stdout=fp )
            return subprocess.call( **kwargs )
    else:
        return subprocess.call( **kwargs )



def do_parallel( func, files, nworkers=None, run=True ):

    # !important - allows one to abort via CTRL-C
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    multiprocessing.log_to_stderr( logging.ERROR )
    
    if not nworkers:
        nworkers = multiprocessing.cpu_count()
    
    p = multiprocessing.Pool( nworkers )

    data = p.map( functools.partial(func, run=run), files )
    p.close()
    p.join()

    return data
