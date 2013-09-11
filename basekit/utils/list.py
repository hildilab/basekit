import collections
import datetime

from basekit.utils.tool import JsonBackend




_ListRecord = collections.namedtuple( "_ListRecord", [
    "name", "source", "retrieved", "version", "list"
])

class ListRecord( _ListRecord ):
    def __init__( self, name, source, retrieved, version, lst ):
        super( ListRecord, self ).__init__(
            name, source, retrieved, version, 
            list( set( map( str.upper, lst ) ) )
        )

class ListIO( JsonBackend ):
    def __init__( self, file_name ):
        super( JsonBackend, self ).__init__( file_name, ListRecord )
    def write( self, records ):
        super( ListIO, self ).write( [ records ] )
    def read( self ):
        return super( ListIO, self ).read()[0]


def list_compare( current_record, compare_record ):
    cur = current_record
    cur_f = frozenset( cur.list )
    cpr = compare_record
    cpr_f = frozenset( cpr.list )
    today = str( datetime.datetime.now().strftime("%Y-%m-%d") )
    new_record = ListRecord(
        "compare", 
        "%s_%s minus %s_%s" % ( 
            cur.name, cur.retrieved, cpr.name, cpr.retrieved 
        ),
        today, today,
        list( cur_f.difference( cpr_f ) )
    )
    old_record = ListRecord(
        "compare", 
        "%s_%s minus %s_%s" % ( 
            cpr.name, cpr.retrieved, cur.name, cur.retrieved
        ),
        today, today,
        list( cpr_f.difference( cur_f ) )
    )
    return new_record, old_record

