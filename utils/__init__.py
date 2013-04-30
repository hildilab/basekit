from __future__ import division
from __future__ import with_statement




def try_int(s, default=None):
    "Convert to integer if possible."
    try: return int(s)
    except: return default or s

def try_float(s, default=None):
    "Convert to float if possible."
    try: return float(s)
    except: return default or s


def get_index(seq, index, default=None):
    if not hasattr(seq, "__getitem__"):
        return default
    try:
        return seq[index]
    except IndexError:
        return default

def boolean(string):
    string = string.lower()
    if string in ['0', 'f', 'false', 'no', 'off']:
        return False
    elif string in ['1', 't', 'true', 'yes', 'on']:
        return True
    else:
        raise ValueError()


def listify( item ):
    """
    Makes a single item a single item list, or returns a list if passed a
    list. Passing a None returns an empty list.
    
    >>> listify( 'a' )
    ['a']
    """
    if not item:
        return []
    elif isinstance( item, list ):
        return item
    else:
        return [ item ]


