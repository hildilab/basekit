from __future__ import division
from __future__ import with_statement


import collections
import itertools
import contextlib
import os
import re



def memoize(f):
    cache = {}
    def memf(*x):
        if x not in cache:
            cache[x] = f(*x)
        return cache[x]
    return memf

def try_int(s, default=None):
    "Convert to integer if possible."
    try: return int(s)
    except: return default if default!=None else s

def try_float(s, default=None):
    "Convert to float if possible."
    try: return float(s)
    except: return default if default!=None else s


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


def iter_overlap( iterator, n=None ):
    if n:
        iterator = itertools.chain( 
            [None]*n, iterator, [None]*n
        )
    return iterator

def iter_window( iterator, n=2, boundary_overlap=None ):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    iterator = iter_overlap( iterator, boundary_overlap )
    result = tuple( itertools.islice( iterator, n ) )
    if len(result)==n:
        yield result    
    for elem in iterator:
        result = result[1:] + (elem,)
        yield result

def iter_stride( iterator, n=2, boundary_overlap=None ):
    "Returns non-overlapping windows (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (sn,s[n+1],...,s[n+n-1]), ...                "
    iterator = iter_overlap( iterator, boundary_overlap )
    while True:
        result = tuple( itertools.islice( iterator, n ) )
        if len(result)!=n:
            return
        yield result

def iter_consume( iterator, n ):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque( iterator, maxlen=0 )
    else:
        # advance to the empty slice starting at position n
        next( itertools.islice( iterator, n, n ), None )


def dir_walker( directory, pattern ):
    for root, dirs, files in os.walk(directory):
        for name in files:
            fpath = os.path.join(root, name)
            m = re.match( pattern, name )
            if( m ):
                yield (m, fpath)


@contextlib.contextmanager 
def working_directory(directory): 
    original_directory = os.getcwd()
    try: 
        os.chdir(directory) 
        yield directory 
    finally: 
        os.chdir(original_directory) 


def copy_dict( dict, **kwargs ):
    dict2 = dict.copy()
    dict2.update( kwargs )
    return dict2
