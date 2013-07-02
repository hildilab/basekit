from __future__ import division

import re
import os
import json
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
    try_int, get_index, boolean, iter_window, memoize, memoize_m,
    copy_dict, dir_walker, DefaultOrderedDict, Bunch, wrap
)
from utils.tool import (
    _, _dir_init, PyTool, DbTool, RecordsMixin, ParallelMixin,
    ProviMixin, SqliteBackend
)

import utils.numpdb as numpdb

from rama import rama_plot
from msa import Muscle


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
    @memoize_m
    def _numa( self, ss ):
        idx_beg = self.npdb.index( 
            chain=ss.chain1, resno=ss.resno1, first=True
        )
        idx_end = self.npdb.index( 
            chain=ss.chain2, resno=ss.resno2, last=True
        )
        return self.npdb.slice( idx_beg, idx_end+1 ).copy( 
            atomname="CA" 
        )
    @memoize_m
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
        _( "angle_next_prev", type="float", default=None, nargs=2 ),
        _( "compare_pdb", type="file", ext="pdb", default=None ),
        _( "compare_sele", type="sele", default=None ),
    ]
    out = [
        _( "query_file", file="query.sql" ),
        _( "provi_dir", dir="provi" ),
        _( "elements_dir", dir="elements" ),
        _( "elements_provi_file", file="elements.provi" ),
        _( "fasta_file", file="elements.fasta" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "sstruc.provi"
    RecordsClass = SstrucDbRecord
    def _init( self, *args, **kwargs ):
        self._init_records( None, **kwargs )
        self.window = range( -2, 3 )
        self.compare_element = None
        self.elements = []
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
        if self.angle_next_prev:
            where.append( 
                "angle_next_prev BETWEEN %i AND %i" % tuple( 
                    self.angle_next_prev 
                )
            )
        where = "\n\tAND ".join( where )
        print "Record count %i" % db.query( where=where, count=True )
        if not self.count:
            self.records = db.query( where=where, limit=self.limit )
        with open( self.query_file, "w" ) as fp:
            fp.write( db.q )
        self.write()
    def _post_exec( self ):
        if self.count:
            return
        phi = DefaultOrderedDict(list)
        psi = DefaultOrderedDict(list)
        
        if self.compare_pdb and self.compare_sele:
            s = self.compare_sele
            b = Bunch(
                chain1=s["chain"], resno1=s["resno"][0],
                chain2=s["chain"], resno2=s["resno"][1],
                pdb_id="COMP", no=0
            )
            n = numpdb.NumPdb( self.compare_pdb )
            self.elements = [( b, n )]

        it = itertools.groupby( 
            self.records, operator.attrgetter('pdb_id') 
        )
        sf = SstrucFinderRefine( 
            [ list(r) for pdb_id, r in it ],
            pdb_archive=self.pdb_archive,
            window=self.window,
            parallel="data"
        )
        self.elements += sf.elements

        # rama_plot( 
        #     zip( phi.values(), psi.values() ),
        #     titles = map( str, window )
        # )
        self._superpose_elements()
        self._msa_elements()
    
    def _superpose_elements( self ):
        if not self.elements:
            return
        provi_elements = []
        r1, numa1 = self.elements[0]
        sele1 = {
            "chain": r1.chain2, "resno": [ r1.resno2-3, r1.resno2 ]
        }
        element_file = os.path.join(
            self.elements_dir, "%s_%i.pdb" % ( r1.pdb_id, r1.no )
        )
        numa1.write( element_file )
        provi_elements.append({ 
            "filename": self.relpath( element_file )
        })
        for r, numa in self.elements[1:]:
            sele = {
                "chain": r.chain2, "resno": [ r.resno2-3, r.resno2 ]
            }
            numpdb.superpose(
                numa, numa1, sele, sele1, align=False
            )
            element_file = os.path.join(
                self.elements_dir, "%s_%i.pdb" % ( r.pdb_id, r.no )
            )
            numa.write( element_file )
            provi_elements.append({ 
                "filename": self.relpath( element_file ),
                "params": { 
                    "load_as": "append",
                    "script": ( "color "
                        "{ chain='%s' and resno>=%i and resno<=%i } "
                        "tomato; " % (r.chain1, r.resno1, r.resno2) )
                }
            })
        with open( self.elements_provi_file, "w" ) as fp:
            json.dump( provi_elements, fp, indent=4 )
    def _msa_elements( self ):
        with open( self.fasta_file, "w" ) as fp:
            for r, numa in self.elements:
                fp.write( ">%s_%i\n%s\n" % (
                    r.pdb_id, r.no, wrap( numa.sequence(), width=80 )
                ))
        Muscle( self.fasta_file )
    




class SstrucFinderRefine( PyTool, ParallelMixin, ProviMixin ):
    args = [
        _( "sstruc_records", type="list" ),
        _( "pdb_archive", type="dir", default="" ),
        _( "window", type="list", default=range(-2, 3) ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "sstruc.provi"
    def _init( self, *args, **kwargs ):
        self._init_parallel( self.sstruc_records, **kwargs )
    def _parallel_results( self, tool_list ):
        self.elements = []
        for t in tool_list:
            self.elements += t.elements
    def func( self ):
        self.sstruc_refined = []
        self.elements = []
        self.phi = []
        self.psi = []
        self.pdb_id = self.sstruc_records[0].pdb_id
        self.pdb_file = self._pdb_file( self.pdb_id )
        self.npdb = numpdb.NumPdb( self.pdb_file )
        for r in self.sstruc_records:
            try:
                C3_O = self.npdb.copy( 
                    chain=r.chain2, resno=r.resno2-2, atomname="O" )
                C2_O = self.npdb.copy( 
                    chain=r.chain2, resno=r.resno2-1, atomname="O" )
                Cd_N = self.npdb.copy( 
                    chain=r.chain2, resno=r.resno2+2, atomname="N" )
                Cdd_N = self.npdb.copy( 
                    chain=r.chain2, resno=r.resno2+3, atomname="N" )
                Cddd_N = self.npdb.copy( 
                    chain=r.chain2, resno=r.resno2+4, atomname="N" )
                if numpdb.numdist( C2_O, Cdd_N ) <= 3.9 or \
                        numpdb.numdist( C3_O, Cd_N ) <= 3.9:
                    self.sstruc_refined.append( r )
            except Exception as e:
                print e
        if self.sstruc_refined:
            self._make_provi()
            # self._make_phi_psi()
            self._make_element()
    def _pdb_file( self, pdb_id ):
        return os.path.join( 
            self.pdb_archive, pdb_id[1:3], "%s.pdb" % pdb_id
        )
    def _make_element( self ):
        for r in self.sstruc_refined:
            idx_beg = self.npdb.index( 
                chain=r.chain1, resno=r.resno1-r.prev_dist, first=True
            )
            idx_end = self.npdb.index( 
                chain=r.chain2, resno=r.resno2+r.next_dist, last=True
            )
            numa = self.npdb.slice( idx_beg, idx_end+1 )
            self.elements.append( ( r, numa ) )
    def _make_provi( self ):
        script = []
        for r in self.sstruc_refined:
            p = (
                r.chain1, r.resno1, r.resno2,
                r.chain1, r.resno1, r.resno2,
                r.chain1, r.resno1-r.prev_dist, r.resno2+r.next_dist,
                r.chain1, r.resno1-r.prev_dist-10, r.resno2+r.next_dist+10,
            )
            s = (
                "center "
                    "{ chain='%s' and resno>=%i and resno<=%i }; "
                "color "
                    "{ chain='%s' and resno>=%i and resno<=%i } tomato; "
                "contact "
                    "{ chain='%s' and resno>=%i and resno<=%i } "
                    "{ not(chain='%s' and resno>=%i and resno<=%i) "
                        "and not water } "
                    "vdw 120%% full;" % p
            )
            script.append( s )
        self._make_provi_file(
            prefix="%s_" % self.pdb_id,
            pdb_file=self.relpath( self.pdb_file ),
            script=" ".join( script )
        )
    def _make_phi_psi( self ):
        for r in self.sstruc_refined:
            for i in self.window:
                numa = self.npdb.copy( chain=r.chain2, resno=r.resno2+i )
                try:
                    self.phi[i].append( numa['phi'][0] )
                    self.psi[i].append( numa['psi'][0] )
                except Exception as e:
                    print i, e
