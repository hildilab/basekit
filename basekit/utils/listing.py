import string
import collections
import datetime
import operator

from basekit import utils
from basekit.utils.tool import JsonBackend


def today():
    return str( datetime.datetime.now().strftime("%Y-%m-%d") )


_ListRecord = collections.namedtuple( "_ListRecord", [
    "name", "source", "retrieved", "version", "length", "list"
])

class ListRecord( _ListRecord ):
    def __new__( _cls, name, source, retrieved, version, lst ):
        if retrieved==None:
            retrieved = today()
        if version==None:
            version = today()
        lst = list( set( map( string.upper, lst ) ) )
        lst.sort()
        return _ListRecord.__new__(
            _cls, name, source, retrieved, version, len(lst), lst
        )

class ListIO( JsonBackend ):
    name = "list"
    def __init__( self, file_name ):
        super( JsonBackend, self ).__init__( file_name, ListRecord )
    def write( self, records ):
        with open( utils.path.mod( self.file_name, ext="txt" ), "w" ) as fp:
            fp.write( "\n".join( records.list ) )
        super( ListIO, self ).write( [ records ] )
    def read( self ):
        return super( ListIO, self ).read()[0]


def list_compare( current_record, compare_record, name=None ):
    cur = current_record
    cur_f = frozenset( cur.list )
    cpr = compare_record
    cpr_f = frozenset( cpr.list )
    if name==None:
        name = "compare"
    new_record = ListRecord(
        name, 
        "%s_%s minus %s_%s" % ( 
            cur.name, cur.retrieved, cpr.name, cpr.retrieved 
        ),
        today(), today(),
        list( cur_f.difference( cpr_f ) )
    )
    old_record = ListRecord(
        name, 
        "%s_%s minus %s_%s" % ( 
            cpr.name, cpr.retrieved, cur.name, cur.retrieved
        ),
        today(), today(),
        list( cpr_f.difference( cur_f ) )
    )
    return new_record, old_record


def list_join( *rec, **kwargs ):
    name = kwargs.get( "name", None )
    if name==None:
        name = "join"
    joined_record = ListRecord(
        name, 
        " plus ".join([ "%s_%s" % ( r.name, r.retrieved ) for r in rec ]),
        today(), today(),
        reduce( operator.add, [ r.list for r in rec ], [] )
    )
    return joined_record


def merge_dic_list(obj_1, obj_2):
    #merges deep dicts with lists in their leafs
    if type(obj_1) == dict and type(obj_2) == dict:
        result = {}
        for key, value in obj_1.iteritems():
            if key not in obj_2:
                result[key] = value
            else:
                result[key] = merge(value, obj_2[key])
        for key, value in obj_2.iteritems():
            if key not in obj_1:
                result[key] = value
        return result
    if type(obj_1) == list and type(obj_2) == list:
        return obj_1 + obj_2
    return obj_2