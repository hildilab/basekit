from __future__ import division
from __future__ import with_statement

import collections
from collections import Hashable
import itertools
import contextlib
import functools
import copy
import os
import re

import job
import path


def class_wraps(cls):
    """Update a wrapper class `cls` to look like the wrapped."""
    class Wrapper(cls):
        """New wrapper that will extend the wrapper `cls` to make it look like `wrapped`.

        wrapped: Original function or class that is beign decorated.
        assigned: A list of attribute to assign to the the wrapper, by default they are:
             ['__doc__', '__name__', '__module__', '__annotations__'].

        """
        def __init__(self, wrapped, assigned=functools.WRAPPER_ASSIGNMENTS):
            print wrapped
            self.__wrapped = wrapped
            for attr in assigned:
                setattr(self, attr, getattr(wrapped, attr))
            super(Wrapper, self).__init__(wrapped)
        def __repr__(self):
            return repr(self.__wrapped)
    return Wrapper


@functools.wraps
def memoize(f):
    cache = {}
    def memf( *args, **kwargs ):
        x = tuple((
            tuple([ 
                v if isinstance(v, Hashable) else tuple(v)
                for v in args
            ]), 
            tuple([ 
                ( k, v if isinstance(v, Hashable) else tuple(v) ) 
                for k, v in kwargs.items()
            ])
        ))
        if x not in cache:
            cache[x] = f(*args, **kwargs)
        return cache[x]
    return memf


@functools.wraps
class memoize_m(object):
    """
    From: http://code.activestate.com/recipes/577452/

    cache the return value of a method
    
    This class is meant to be used as a decorator of methods. The return value
    from a given method invocation will be cached on the instance whose method
    was invoked. All arguments more passed to a method decorated with memoize must
    be hashable.
    
    If a memoized method is invoked directly on its class the result will not
    be cached. Instead the method will be invoked like a static method:
    class Obj(object):
        @memoize
        def add_to(self, arg):
            return self + arg
    Obj.add_to(1) # not enough arguments
    Obj.add_to(1, 2) # returns 3, result is not cached
    """
    def __init__(self, func):
        self.func = func
    def __get__(self, obj, objtype=None):
        if obj is None:
            return self.func
        return functools.partial(self, obj)
    def __call__(self, *args, **kwargs):
        obj = args[0]
        try:
            cache = obj.__cache
        except AttributeError:
            cache = obj.__cache = {}
        key = (
            self.func,
            tuple([ 
                v if isinstance(v, Hashable) else tuple(v)
                for v in args[1:]
            ]),
            tuple([ 
                ( k, v if isinstance(v, Hashable) else tuple(v) ) 
                for k, v in kwargs.items()
            ])
        )
        try:
            res = cache[key]
        except KeyError:
            res = cache[key] = self.func(*args, **kwargs)
        return res


@functools.wraps
def memoize1(f):
    """ Memoization decorator for a function taking a single argument """
    class memodict(dict):
        def __missing__(self, key):
            ret = self[key] = f(key)
            return ret 
    return memodict().__getitem__


def lazy_property(fn):
    attr_name = '_lazy_' + fn.__name__
    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazy_property

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

def wrap(text, width):
    """
    http://code.activestate.com/recipes/148061-one-liner-word-wrap-function/
    A word-wrap function that preserves existing line breaks
    and most spaces in the text. Expects that existing line
    breaks are posix newlines (\n).
    """
    return reduce(lambda line, word, width=width: '%s%s%s' %
                  (line,
                   ' \n'[(len(line)-line.rfind('\n')-1
                         + len(word.split('\n',1)[0]
                              ) >= width)],
                   word),
                  text.split(' ')
                 )

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

def flatten( lists ):
    return list( itertools.chain.from_iterable( lists ) )

def iter_overlap( iterator, n=None ):
    iterator = iter(iterator)
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
    iterator = iter(iterator)
    # Use functions that consume iterators at C speed.
    if n is None:
        # feed the entire iterator into a zero-length deque
        collections.deque( iterator, maxlen=0 )
    else:
        # advance to the empty slice starting at position n
        next( itertools.islice( iterator, n, n ), None )

def dir_walker( directory, pattern ):
    for root, dirs, files in os.walk(directory):
        files.sort()
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



class Bunch( object ):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


class DefaultOrderedDict( collections.OrderedDict ):
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory!=None and not callable(default_factory)):
            raise TypeError('first argument must be callable')
        collections.OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory
    def __getitem__(self, key):
        try:
            return collections.OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value
    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()
    def copy(self):
        return self.__copy__()
    def __copy__(self):
        return type(self)(self.default_factory, self)
    def __deepcopy__(self, memo):
        return type(self)(
            self.default_factory, copy.deepcopy(self.items())
        )
    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (
            self.default_factory, collections.OrderedDict.__repr__(self)
        )

