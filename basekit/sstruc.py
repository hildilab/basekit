#! /usr/bin/env python

from __future__ import division

import re
import os
import argparse
import operator
import sqlite3
import functools
import itertools
import collections
import string
import textwrap

import numpy as np

import utils.math
import utils.path
from utils import (
    try_int, get_index, boolean, iter_window, memoize, 
    copy_dict, dir_walker
)
from utils.timer import Timer
from utils.db import get_pdb_files, create_table
from utils.tool import (
    _, _dir_init, PyTool, DbTool, RecordsMixin, ParallelMixin,
    ProviMixin, SqliteBackend
)

import utils.numpdb as numpdb

from rama import rama_plot

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "sstruc" )

SstrucDbRecord = collections.namedtuple( 'SstrucDbRecord', [
    'pdb_id', 'type', 'subtype', 'length',
    'chain1', 'resno1', 'chain2', 'resno2',
    'prev_type', 'prev_dist', 'next_type', 'next_dist',
    'hbond_cur_prev', 'hbond_next_cur', 'hbond_next_prev',
    'angle_cur_prev', 'angle_next_cur', 'angle_next_prev',
    'no'
])


# TODO length calculations do not take altloc into account
# (which will be handled in numpdb) but also no missing residues
# (maybe add fake resiudes to numpdb for that purpose???)
class BuildSstrucDbRecords( object ):
    def __init__( self, pdb_file, pdb_id=None ):
        self.npdb = numpdb.NumPdb( pdb_file, features={ "phi_psi": False } )
        self.pdb_id = pdb_id or utils.path.stem( pdb_file )
        # for a in self.npdb._atoms: print a
        # for a in self.npdb._iter_resno(): print a._atoms
    def _sheet_hbond( self, x, y ):
        if not x or not y:
            return None
        if x.type!=numpdb.SHEET or y.type!=numpdb.SHEET:
            return None
        # TODO check if in idx interval, but probably ok as is
        # because there is only a resno given, no chain, altloc, insertion
        return y.resno1 <= x.hbond <= y.resno2
    @memoize
    def _numa( self, ss ):
        idx_beg = self.npdb.index( 
            chain=ss.chain1, resno=ss.resno1, interval=True
        )
        idx_end = self.npdb.index( 
            chain=ss.chain2, resno=ss.resno2, interval=True
        )
        return self.npdb.slice( idx_beg[0], idx_end[-1]+1 ).copy( 
            atomname="CA" 
        )
    @memoize
    def _axis( self, ss ):
        numa = self._numa( ss )
        if numa.length<3:
            return None
        beg, end = numa.axis()
        return end-beg
    def _sstruc_angle( self, x, y ):
        if not x or not y:
            return None
        ax = self._axis( x )
        ay = self._axis( y )
        if ax==None or ay==None:
            return None
        else:
            return round( utils.math.angle( ax, ay ), 3 )
    def _length( self, ss ):
        numa = self._numa( ss )
        return numa.length
    def _distance( self, ss1, ss2 ):
        if not ss1 or not ss2:
            return None
        if ss1.chain1!=ss2.chain1:
            return None
        idx_beg = self.npdb.index( 
            chain=ss1.chain2, resno=ss1.resno2, interval=True
        )
        idx_end = self.npdb.index( 
            chain=ss2.chain1, resno=ss2.resno1, interval=True
        )
        if idx_beg==idx_end:
            return -1
        if idx_beg[-1]+1==idx_end[0]:
            return 0
        # overlap
        if idx_beg>idx_end:
            numa = self.npdb.slice( idx_end[0], idx_beg[-1]+1 ).copy( 
                atomname="CA" 
            )
            return numa.length * -1
        else:
            numa = self.npdb.slice( idx_beg[-1]+1, idx_end[0] ).copy( 
                atomname="CA" 
            )
            return numa.length
    def _make_record( self, i, cur, prev, next ):
        return SstrucDbRecord(
            self.pdb_id, cur.type, cur.subtype,
            self._length( cur ),
            cur.chain1, cur.resno1,
            cur.chain2, cur.resno2,
            prev and prev.type,
            self._distance( prev, cur ),
            next and next.type,
            self._distance( cur, next ),
            self._sheet_hbond( cur, prev ),
            self._sheet_hbond( next, cur ),
            self._sheet_hbond( next, prev ),
            self._sstruc_angle( cur, prev ),
            self._sstruc_angle( next, cur ),
            self._sstruc_angle( next, prev ),
            i
        )
    def get( self ):
        records = []
        it = iter_window( self.npdb._sstruc, 3, boundary_overlap=1 )
        for i, d in enumerate( it, start=1 ):
            prev, cur, next = d
            records.append( self._make_record( i, cur, prev, next ) )
        return records



class Sstruc( PyTool, RecordsMixin, ParallelMixin ):
    args = [
        _( "pdb_input", type="file", ext="pdb" ),
        _( "pdb_id", type="text", default=None )
    ]
    RecordsClass = SstrucDbRecord
    def _init( self, *args, **kwargs ):
        self._init_records( self.pdb_input, **kwargs )
        self._init_parallel( self.pdb_input, **kwargs )
    def func( self ):
        self.records = BuildSstrucDbRecords( self.pdb_input ).get()
        self.write()


class SstrucInfo( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "sele", type="sele" )
    ]
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        npdb = numpdb.NumPdb( self.pdb_file )
        self.sele[ "atomname" ] = "CA"
        numa = npdb.copy( **self.sele )
        sstruc_order = []
        for k, g in itertools.groupby( numa['sstruc'] ):
            n = len( list(g) )
            sstruc_order.append( ( k, n ) )
            print k, n


class SstrucFinder( PyTool, RecordsMixin, ProviMixin ):
    args = [
        _( "sstruc_db", type="file", ext="sqlite" ),
        _( "pdb_archive", type="dir", default=None ),
        _( "count", type="checkbox", default=False ),
        _( "limit", type="int", default=None ),
        _( "type", type="str", default="H", help="H or E" ),
        _( "subtype", type="int", default=None, help="H: ; E: " ),
        _( "length", type="int", default=None, nargs=2, help="" ),
        _( "prev_type", type="str", default=None, help="H or E" ),
        _( "prev_dist", type="int", default=None, nargs=2, help="" ),
        _( "next_type", type="str", default=None, help="H or E" ),
        _( "next_dist", type="int", default=None, nargs=2, help="" ),
        _( "next_dist", type="int", default=None, nargs=2, help="" ),
        _( "hbond_next_prev", type="checkbox", default=False ),
    ]
    out = [
        _( "query_file", file="query.sql" ),
        _( "provi_dir", dir="provi" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "sstruc.provi"
    RecordsClass = SstrucDbRecord
    def _init( self, *args, **kwargs ):
        self._init_records( None, **kwargs )
    def func( self ):
        db = SqliteBackend( self.sstruc_db, SstrucDbRecord )
        where = []
        types = { "H": 1, "E": 2 }
        if self.type:
            where.append( "type=%i" % types[ self.type ] )
        if self.subtype:
            where.append( "subtype=%i" % self.subtype )
        if self.length:
            where.append( 
                "length BETWEEN %i AND %i" % tuple(self.length) 
            )
        if self.prev_type:
            where.append( "prev_type=%i" % types[ self.prev_type ] )
        if self.next_type:
            where.append( "next_type=%i" % types[ self.next_type ] )
        if self.prev_dist:
            where.append( 
                "prev_dist BETWEEN %i AND %i" % tuple( self.prev_dist )
            )
        if self.next_dist:
            where.append( 
                "next_dist BETWEEN %i AND %i" % tuple( self.next_dist ) 
            )
        if self.hbond_next_prev:
            where.append( "hbond_next_prev=1" )
        where = "\n\tAND ".join( where )
        if self.count:
            print "Count %i" % db.query( where=where, count=True )
        else:
            self.records = db.query( where=where, limit=self.limit )
        with open( self.query_file, "w" ) as fp:
            fp.write( db.q )
        self.write()
    def _post_exec( self ):
        if self.count:
            return
        pdb_groups = []
        phi = []
        psi = []
        it = itertools.groupby( 
            self.records, operator.attrgetter('pdb_id') 
        )
        for pdb_id, records in it:
            rec = list(records)
            self._make_provi( pdb_id, rec )
            _phi, _psi = self._get_phi_psi( pdb_id, rec )
            phi += _phi
            psi += _psi
        self._make_rama( phi, psi )
    def _pdb_file( self, pdb_id ):
        return os.path.join( 
            self.pdb_archive, pdb_id[1:3], "%s.pdb" % pdb_id
        )
    def _make_provi( self, pdb_id, records ):
        pdb_file = self._pdb_file( pdb_id )
        script = []
        for r in records:
            s = "color { chain='%s' and %s-%s } tomato;" % ( 
                r.chain1, r.resno1, r.resno2 
            )
            script.append( s )
        self._make_provi_file(
            output_dir=self.provi_dir,
            prefix="%s_" % pdb_id,
            pdb_file=os.path.relpath( pdb_file, self.provi_dir ),
            script=" ".join( script )
        )
    def _get_phi_psi( self, pdb_id, records ):
        npdb = numpdb.NumPdb( self._pdb_file( pdb_id ) )
        phi = []
        psi = []
        for r in records:
            numa = npdb.copy( chain=r.chain2, resno=r.resno2 )
            phi.append( numa['phi'][0] )
            psi.append( numa['psi'][0] )
        return phi, psi
    def _make_rama( self, phi, psi ):
        rama_plot( ( phi, psi ) )


        




def sstruc_test( pdb_file ):
    pass

class SstrucTest( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    no_output = True
    def func( self ):
        sstruc_test( self.pdb_file )





def sstruc2jmol( sstruc ):
    ret = ""
    for i, ss in enumerate( sstruc ):
        ret += "draw ID 'v%i' vector {%s} {%s};" % (
            i,
            "%0.2f %0.2f %0.2f" % tuple(ss[8]),
            "%0.2f %0.2f %0.2f" % tuple(ss[7])
        )
    return ret




def create_provi( output_dir, template, row ):
    pdb_path = row[0]
    sele = "{ chain='%s' and %s-%s }" % ( row[3], row[4], row[6] )
    pdb_fdir, pdb_fname = os.path.split( pdb_path )
    fname = '%s_%s%s-%s.provi' % ( pdb_fname[0:4], row[3], row[4], row[6] )
    fpath = os.path.join( output_dir, fname )
    npdb = numpdb.NumPdb( pdb_path )
    d = {
        "pdb_file": pdb_path,
        "script": "select %s; color {selected} tomato; center {selected};%s" % ( 
            sele, 
            numpdb.sstruc2jmol( 
                filter( lambda x: x[10] in [row[18]-1, row[18], row[18]+1], npdb.sstruc )
            ) 
        )
    }
    create_json( d, fpath, template )


