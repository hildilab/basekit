from __future__ import division
from __future__ import with_statement

import contextlib
import os
import re



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


def iter_window(iterator, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    #iterator = iter(seq)
    result = tuple(itertools.islice(iterator, n))
    if len(result) == n:
        yield result    
    for elem in iterator:
        result = result[1:] + (elem,)
        yield result


def iter_consume(iterator, n):
    "Advance the iterator n-steps ahead. If n is none, consume entirely."
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque(iterator, maxlen=0)
    else:
        # advance to the empty slice starting at position n
        next(itertools.islice(iterator, n, n), None)


def dir_walker( dir, pattern ):
    for root, dirs, files in os.walk(top):
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

