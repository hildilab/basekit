

import os
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
        stem( file_name ) + suffix + ext( filename )
    )

def change_ext( file_name, extension="" ):
	return join( 
        dirname( file_name ),
        stem( file_name ) + "." + extension
    )






