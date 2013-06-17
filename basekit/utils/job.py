
import multiprocessing
import os
import sys
import subprocess
import fcntl
import select
import time
import signal
import logging

try:
    from subprocess import DEVNULL # python 3.3
except ImportError:
    DEVNULL = open( os.devnull, 'wb' )



def run_command( cmd, cwd=".", log=None, verbose=False ):
    kwargs = {
        "args": map( str, cmd ),
        "cwd": cwd,
        "stdout": DEVNULL,
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



def run_parallel( cmd_list, nworkers=None, verbose=False ):
    """UNTESTED"""

    # !important - allows one to abort via CTRL-C
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    multiprocessing.log_to_stderr( logging.ERROR )
    
    if not nworkers: nworkers = multiprocessing.cpu_count()
    p = multiprocessing.Pool( nworkers, maxtasksperchild=50 )

    data = []

    for cmd in cmd_list:
        p.apply_async( 
            cmd, 
            kwds={ "verbose": verbose },
            callback=data.append
        )
    
    p.close()
    p.join()

    return data
