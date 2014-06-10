

import os
import errno
from os.path import basename, splitext, join, dirname


def ext( file_name ):
    return splitext( file_name )[1]


def stem( file_name ):
    return splitext( basename( file_name ) )[0]


def add_prefix( file_name, prefix="" ):
    return join(
        dirname( file_name ),
        prefix + basename( file_name )
    )


def add_suffix( file_name, suffix="" ):
    return join(
        dirname( file_name ),
        stem( file_name ) + suffix + ext( file_name )
    )


def change_ext( file_name, extension="" ):
    return join(
        dirname( file_name ),
        stem( file_name ) + "." + extension
    )


def change_name( file_name, name="" ):
    return join(
        dirname( file_name ),
        name + ext( file_name )
    )


def mod( file_name, ext=None, name=None, prefix=None, suffix=None ):
    if ext:
        file_name = change_ext( file_name, ext )
    if name:
        file_name = change_name( file_name, name )
    if prefix:
        file_name = add_prefix( file_name, prefix )
    if suffix:
        file_name = add_suffix( file_name, suffix )
    return file_name


def remove( filename ):
    """silently remove a file"""
    try:
        os.remove( filename )
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise e


def which( program ):
    # from http://stackoverflow.com/a/377028
    def is_exe( fpath ):
        return os.path.isfile( fpath ) and os.access( fpath, os.X_OK )
    fpath, fname = os.path.split( program )
    if fpath:
        if is_exe( program ):
            return program
    else:
        for path in os.environ["PATH"].split( os.pathsep ):
            path = path.strip( '"' )
            exe_file = os.path.join( path, program )
            if is_exe( exe_file ):
                return exe_file
    return None
