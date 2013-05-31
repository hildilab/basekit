

import os
import string



def stem( file_name ):
    return os.path.splitext( os.path.split( file_name )[-1] )[0]


def tmpl( template, file_name, directory=None ):
    values = { 
    	"stem": stem( file_name )
	}
    new_name = string.Template( template ).substitute( **values )
    if directory:
        new_name = os.path.join( directory, new_name )
    return new_name


def new_ext( file_name, new_ext, directory=None ):
    return tmpl( "${stem}." + new_ext )






