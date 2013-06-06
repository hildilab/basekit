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
import multiprocessing
import signal
import logging


import numpy as np

import utils.math
import utils.path
from utils import (
    try_int, get_index, boolean, iter_window, memoize, 
    copy_dict, dir_walker
)
from utils.timer import Timer
from utils.db import get_pdb_files, create_table
from utils.tool import PyTool, DbTool, RecordsMixin, ParallelMixin

import utils.numpdb as numpdb


logging.basicConfig()
LOG = logging.getLogger('sstruc')
# LOG.setLevel( logging.ERROR )
LOG.setLevel( logging.WARNING )
# LOG.setLevel( logging.DEBUG )


# http://docs.python.org/2/library/collections.html#collections.namedtuple


SstrucDbRecord = collections.namedtuple( 'SstrucDbRecord', [
    'pdb_id', 'type', 'subtype', 'length',
    'chain1', 'resno1', 'chain2', 'resno2',
    'prev_type', 'prev_dist', 'next_type', 'next_dist',
    'hbond_cur_prev', 'hbond_next_cur', 'hbond_next_prev',
    'angle_cur_prev', 'angle_next_cur', 'angle_next_prev',
    'no'
])


class BuildSstrucDbRecords( object ):
    def __init__( self, pdb_file, pdb_id=None ):
        self.npdb = numpdb.NumPdb( pdb_file, features={ "phi_psi": True } )
        self.pdb_id = pdb_id
        # for a in self.npdb._atoms: print a
        # for a in self.npdb._iter_resno(): print a._atoms
    def _sheet_hbond( self, x, y ):
        if not x or not y:
            return None
        if x.type!=numpdb.SHEET or y.type!=numpdb.SHEET:
            return None
        return y.resno1 <= x.hbond <= y.resno2
    @memoize
    def _axis( self, ss ):
        idx_beg = self.npdb.index(
            chain=ss.chain1, resno=ss.resno1, first=True 
        )
        idx_end = self.npdb.index( 
            chain=ss.chain2, resno=ss.resno2, last=True
        )
        if idx_beg==idx_end:
            # no axis for a point
            return None
        else:
            if idx_beg>idx_end:
                idx_beg, idx_end = idx_end, idx_beg
            beg, end = self.npdb.slice( idx_beg, idx_end+1 ).axis(
                atomname=numpdb.ATOMS["backbone"]
            )
            return end-beg
    def _sstruc_angle( self, x, y ):
        if not x or not y:
            return None
        ax = self._axis( x )
        ay = self._axis( y )
        if ax==None or ay==None:
            return None
        else:
            return utils.math.angle( ax, ay )
    def _make_record( self, i, cur, prev, next ):
        return SstrucDbRecord(
            self.pdb_id, cur.type, cur.subtype,
            cur.resno2 - cur.resno1 + 1,
            cur.chain1, cur.resno1,
            cur.chain2, cur.resno2,
            prev and prev.type,
            prev and cur.resno1 - prev.resno2,
            next and next.type,
            next and next.resno1 - cur.resno2,
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
        { "name": "pdb_input", "type": "file", "ext": "pdb" },
        { "name": "pdb_id", "type": "text", "default_value": None }
    ]
    RecordsClass = SstrucDbRecord
    def _init( self, pdb_input, pdb_id=None, **kwargs ):
        self.pdb_input = self.abspath( pdb_input )
        self.pdb_id = pdb_id
        self.id = pdb_id
        self.output_files = []
        self._init_records( utils.path.stem( pdb_input ), **kwargs )
        self._init_parallel( self.pdb_input, **kwargs )
        #print self.records[0] if self.records else None
    def func( self ):
        if self.parallel:
            self._make_tool_list()
            tool_results = self._func_parallel()
            self.records = list(itertools.chain.from_iterable(
                map( operator.attrgetter( "records" ), tool_results )
            ))
            # print self.pdb_input, utils.path.stem( self.pdb_input )
            # print list(self.records)[0] if self.records else None
        else:
            self.records = BuildSstrucDbRecords( 
                self.pdb_input, pdb_id=self.pdb_id 
            ).get()
            # print self.records[0] if self.records else None
        self.write()
        # for r in self.records:
        #     print r
        # with Timer( "write sstruc: %s" % self.output_type ):
        #     self.write()
        # print self.records[0]
        # self.records = None
        # with Timer( "read sstruc: %s" % self.output_type ):
        #     self.read()
        # print self.records[0]



class SstrucDb( DbTool ):
    pass


class SstrucFinder( PyTool ):
    pass


class SstrucPlot( PyTool ):
    pass




def sstruc_test( pdb_file ):
    pass

class SstrucTest( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    no_output = True
    def _init( self, pdb_file, **kwargs ):
        self.pdb_file = self.abspath( pdb_file )
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


def _sstruc_pdb( fpath ):
    return get_sstruc_records( numpdb.NumPdb( fpath ) )


def sstruc_pdb( pdb_files, nworkers=None ):
    return do_parallel( _sstruc_pdb, pdb_files, nworkers )


def create_sstruc_db( db_path, data, overwrite=False ):
    schema = textwrap.dedent("""
        CREATE TABLE sstruc (
            pdb_id text, type int, subtype int, 
            length int, chain1 text, resno1 int, chain2 text, resno2 int,
            prev_type int, prev_dist int, next_type int, next_dist int,
            cur_prev_hbond boolean, cur_next_hbond boolean, prev_next_hbond boolean,
            cur_prev_angle real, cur_next_angle real, prev_next_angle real,
            no int
        )"""
    )
    name = "sstruc"

    all_data = []
    for pdb_path, sstruc_record in data:
        if sstruc_record:
            all_data += [ [pdb_path] + r for r in sstruc_record ]

    return create_table( db_path, schema, name, all_data, overwrite=overwrite )



def query_sstruc_db( db_path, query, func=None ):
    with sqlite3.connect( db_path, isolation_level="EXCLUSIVE" ) as conn:
        c = conn.cursor()
        for row in c.execute( query ):
            if func:
                func( row )
            
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




def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-db_file', help='path to the db', type=str)
    parser.add_argument(
        '-create_db', help='create db', type=boolean)
    parser.add_argument(
        '-overwrite', help='overwrite existing table when creating db', type=boolean)
    parser.add_argument(
        '-pdb_path', help='path to the pdb files and directories', type=str)
    parser.add_argument(
        '-test_sample', help='how many?', type=int)
    parser.add_argument(
        '-query_file', help='path to a sql file', type=str)
    parser.add_argument(
        '-provi_tpl', help='path to a provi template', type=str)
    parser.add_argument(
        '-provi_out', help='path where the provi files should be written to', type=str)
    parser.add_argument(
        '-analyze', help='do some analysis', type=boolean)
    parser.add_argument(
        '-pdb', help='get sstruc record for a pdb file', type=str)
    parser.add_argument(
        '-dssp', help='parse dssp for a pdb file', type=str)

    
    # parse the command line
    args = parser.parse_args()



    if args.pdb_path:
        pdb_files = get_pdb_files( args.pdb_path, pattern=".pdb" )
        if args.test_sample:
            pdb_files = pdb_files[0:args.test_sample]
        print "%s pdb files" % len( pdb_files )


    if args.create_db and args.db_file and args.pdb_path: 
        with Timer("get sstruc pdb -> records"):
            sstruc_records = sstruc_pdb( pdb_files )

        # with Timer("struc records"):
        #     sstruc_records = [ (x[0], get_sstruc_records(x[1])) for x in sstruc_npdb_list ]

        with Timer("create sstruc db"):
            create_sstruc_db( args.db_file, sstruc_records, overwrite=args.overwrite )


    if args.query_file and args.db_file:
        with open( args.query_file, "r" ) as fp:
            query = fp.read()
        if args.provi_tpl and args.provi_out:
            if not os.path.exists( args.provi_out ):
                os.makedirs( args.provi_out )
            fn = functools.partial( create_provi, args.provi_out, args.provi_tpl )
            query_sstruc_db( args.db_file, query, fn )


    if args.analyze and args.pdb_path:
        with Timer("pdb2numpy"):
            for fpath in pdb_files:
                #print fpath
                npdb = numpdb.NumPdb( fpath )
                for i, ss in enumerate( npdb.sstruc() ):
                    # print i, ss[4]-ss[2]+1
                    a = npdb.axis( chain=ss[1], atomname="CA", resno=[ ss[2], ss[4] ] )
                    print "draw ID 'v%i' vector {%s} {%s};" % (
                        i,
                        "%0.2f %0.2f %0.2f" % tuple(a[0]),
                        "%0.2f %0.2f %0.2f" % tuple(a[1]-a[0])
                    )


    if args.pdb:
        npdb = numpdb.NumPdb( args.pdb )
        records = get_sstruc_records( npdb )
        for r in records:
            print r
            print npdb.sequence( resno=[r[3], r[5]+4], chain=r[2], atomname="CA" )


    if args.dssp:
        s = parse_dssp( args.dssp )
        with open( "%s.txt" % args.dssp, "w" ) as fp:
            fp.write( s )


if __name__ == "__main__":
    main()