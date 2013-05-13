#! /usr/bin/env python

from __future__ import division

import re
import os
import argparse
import operator
import sqlite3
import functools

from collections import defaultdict
from string import Template

import numpy as np

from utils import try_int, get_index, boolean
from utils.timer import Timer
from utils.db import get_pdb_files, create_table
from utils.job import do_parallel
from utils.math import vec_angle

import utils.numpdb as numpdb



def _sstruc_pdb( fpath ):
    return get_sstruc_records( numpdb.NumPdb( fpath ) )


def sstruc_pdb( pdb_files, nworkers=None ):
    return do_parallel( _sstruc_pdb, pdb_files, nworkers )

def _sheet_test( x, y ):
    if x[0]!=numpdb.SHEET or y[0]!=numpdb.SHEET:
        return None
    return y[2] <= x[6] <= y[4]


def make_struc_record( cur, prev, next ):
    return [
        cur[0],                 # type
        cur[4] - cur[2] + 1,    # length
        cur[1],                 # chain 1
        cur[2],                 # resno 1
        cur[3],                 # chain 2
        cur[4],                 # resno 2
        cur[5],                 # subtype
        prev and prev[0],           # prev type
        prev and cur[2] - prev[4],  # prev dist
        next and next[0],           # next type
        next and next[2] - cur[4],  # next dist
        prev and _sheet_test( cur, prev ),              # hbond cur prev
        next and _sheet_test( next, cur ),              # hbond next cur
        prev and next and _sheet_test( next, prev ),    # hbond next prev
        prev and vec_angle( cur[7], prev[7] ),              # angle cur prev
        next and vec_angle( next[7], cur[7] ),              # angle next cur
        prev and next and vec_angle( next[7], prev[7] ),    # angle next prev
        cur[10]                 # running number
    ]


def _sstruc_get( x, j, sstruc ):
    if j<0 or j>=len(sstruc): return None
    y = sstruc[j]
    if x[1]==y[1]:
        return y
    else:
        return None

def get_sstruc_records( npdb ):
    ss = npdb.sstruc
    r = []
    for i in range(len(ss)):
        cur = ss[i]
        prev = _sstruc_get( cur, i-1, ss )
        next = _sstruc_get( cur, i+1, ss )
        r.append( make_struc_record( cur, prev, next ) )
    return r




def create_sstruc_db( db_path, data, overwrite=False ):
    schema = "CREATE TABLE sstruc (pdb_path text, type int, length int, chain1 text, resno1 int, chain2 text, resno2 int, subtype int, prev_type int, prev_dist int, next_type int, next_dist int, cur_prev_hbond boolean, cur_next_hbond boolean, prev_next_hbond boolean, cur_prev_angle real, cur_next_angle real, prev_next_angle real, no int)"
    name = "sstruc"

    all_data = []
    for pdb_path, sstruc_record in data:
        if sstruc_record:
            all_data += [ [pdb_path] + r for r in sstruc_record ]

    return create_table( db_path, schema, name, all_data, overwrite=overwrite )


def create_json( values_dict, out, tpl ):
    with open( tpl, "r" ) as fp:
        tpl_str = fp.read()
    with open( out, "w" ) as fp:
        fp.write( Template( tpl_str ).substitute( **values_dict ) )


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


# http://swift.cmbi.ru.nl/gv/dssp/HTML/descrip.html
def parse_dssp( dssp_path ):
    records = [ None ]
    with open( dssp_path, "r" ) as fp:
        residue_section = False
        for line in fp:
            if line.startswith( "  #  RESIDUE AA STRUCTURE" ):
                residue_section = True
            elif line.strip() and residue_section:
                if line[13]=="!":
                    records.append( None )
                else:
                    records.append([
                        int(line[0:5]),                 # no
                        int(line[5:10]),                # resno
                        line[11],                       # chain
                        line[16],                       # ss type
                        int(line[25:29]),               # hbond to resno
                        line[13],                       # aa1
                        line[33],                       # beta sheet label
                    ])
    prev_ss = " "
    prev_chain = " "
    sstruc = [[]]
    for i, r in enumerate( records ):
        if r:
            if r[3]==prev_ss:
                sstruc[-1].append(r)
            else:
                sstruc.append([r])
            prev_chain = r[2]
            prev_ss = r[3]
    sheet_label_count = defaultdict(int)
    for ss in sstruc:
        if ss[0][3] in [ "E" ]:
            sheet_label_count[ ss[0][6] ] += 1
    helix_counter = 0
    helix_list = []
    sheet_counter = 0
    sheet_list = []
    for ss in sstruc:
        if ss[0][3] in [ "H" ]:
            helix_counter += 1
            helix_list.append(
                "HELIX% 5i% 4i %s %s% 5i  %s %s% 5i % 2i                               % 5i    " % (
                    helix_counter, helix_counter, 
                    numpdb.AA3.get( ss[0][5], "XXX" ), ss[0][2], ss[0][1], 
                    numpdb.AA3.get( ss[-1][5], "XXX" ), ss[-1][2], ss[-1][1],
                    1, len(ss)
                )
            )
        if ss[0][3] in [ "E" ]:
            sheet_counter += 1
            e = records[ ss[0][4] ]
            e2 = records[ ss[-1][4] ]
            if not e2 or ss[0][0]<e2[0]:
                hbond = "                                        "
            else:
                hbond = "  O  %s %s% 4i   N  %s %s% 4i           " % (
                    numpdb.AA3.get( ss[-1][5], "XXX" ), ss[-1][2], ss[-1][1],
                    numpdb.AA3.get( e2[5], "XXX" ), e2[2], e2[1],
                )
            sheet_list.append(
                "SHEET% 5i   %s% 2i %s %s% 4i  %s %s% 4i -1%s" % (
                    sheet_counter, ss[0][6], sheet_label_count[ ss[0][6] ],
                    numpdb.AA3.get( ss[0][5], "XXX" ), ss[0][2], ss[0][1], 
                    numpdb.AA3.get( ss[-1][5], "XXX" ), ss[-1][2], ss[-1][1],
                    hbond
                )
            )
    # print "HELIX    1   1 GLU A   33  HIS A   65  1                                  33    "
    # for h in helix_list:
    #     print h
    # print "SHEET    2   A 2 TYR A  10  VAL A  11 -1  O  VAL A  11   N  THR A   4           "
    # for e in sheet_list:
    #     print e
    return "\n".join( helix_list + sheet_list )


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