#! /usr/bin/env python

from __future__ import division

import re
import os
import argparse
import gzip

import operator
import sqlite3
import cPickle

from utils import try_int, get_index, boolean
from utils.timer import Timer

from utils.db import get_pdb_files

from utils.job import _prep_func, do_parallel

HELIX = 1
SHEET = 2







     
@_prep_func
def _unzip_pdb( fpath ):
    fdir, fname = os.path.split( fpath )
    fpdb = os.path.join( fdir, fname[3:7]+".pdb" )
    with open( fpdb, "w" ) as f:
        with gzip.open( fpath, "rb" ) as z:
            f.write( z.read() )
    return True

def unzip_pdb( pdb_files, nworkers=None ):
    return do_parallel( _unzip_pdb, pdb_files, nworkers )


@_prep_func
def _sstruc_pdb( fpath ):
    sstruc = []
    with open( fpath, "r" ) as f:
        for line in f:
            if line.startswith("HELIX"):
                sstruc.append([
                    HELIX,
                    line[19],                   # chain 1
                    try_int( line[21:25] ),     # resno 1
                    line[31],                   # chain 2
                    try_int( line[33:37] ),     # resno 2
                    try_int( line[38:40] ),     # subtype
                ])
            elif line.startswith("SHEET"):
                sstruc.append([
                    SHEET,
                    line[21],                   # chain 1
                    try_int( line[22:26] ),     # resno 1
                    line[32],                   # chain 2
                    try_int( line[33:37] ),     # resno 2
                    try_int( line[38:40] ),     # strand sense (subtype)
                    try_int( line[65:69] ),     # resno hbond prev strand
                ])
            elif line.startswith("ATOM"):
                break
    sstruc.sort( key=operator.itemgetter(1,2) )
    return sstruc


def sstruc_pdb( pdb_files, nworkers=None ):
    return do_parallel( _sstruc_pdb, pdb_files, nworkers )

def _sheet_test( x, y ):
    if x[0]!=SHEET or y[0]!=SHEET:
        return None
    return y[2] <= x[6] <= y[4]


def make_struc_record( cur, prev, next ):
    return [
        cur[0],             # type
        cur[4] - cur[2],    # length
        cur[1],             # chain 1
        cur[2],             # resno 1
        cur[3],             # chain 2
        cur[4],             # resno 2
        cur[5],             # subtype
        prev and prev[0],           # prev type
        prev and cur[2] - prev[4],  # prev dist
        next and next[0],           # next type
        next and next[2] - cur[4],   # next dist
        prev and _sheet_test( cur, prev ),
        next and _sheet_test( next, cur ),
        prev and next and _sheet_test( next, prev ),
    ]

def _sstruc_get( x, j, data ):
    if j<0: return None
    y = get_index( data, j )
    if y and x[1]==y[1]:
        return y
    else:
        return None

def get_sstruc_records( data ):
    r = []
    for i in range(len(data)):
        cur = data[i]
        prev = _sstruc_get( cur, i-1, data )
        next = _sstruc_get( cur, i+1, data )
        r.append( make_struc_record( cur, prev, next ) )
    return r

def create_sstruc_db( data ):
    conn = sqlite3.connect('sstruc.db', isolation_level="EXCLUSIVE")
    c = conn.cursor()

    c.execute('''DROP TABLE sstruc''')

    # create table
    c.execute('''CREATE TABLE sstruc
                (pdb_path text, type int, length int, chain1 text, resno1 int, chain2 text, resno2 int, subtype int, prev_type int, prev_dist int, next_type int, next_dist int, cur_prev_hbond boolean, cur_next_hbond boolean, prev_next_hbond boolean)''')

    # populate table
    all_data = []
    for pdb, sstruc in data:
        all_data += [ [pdb] + r for r in sstruc ]
        
    c.executemany('INSERT INTO sstruc VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', all_data)

    conn.commit()
    conn.close()


# select * from sstruc where type="HELIX" and prev_type="SHEET" and next_type="SHEET"


def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-init', help='initialize db', type=boolean)

    
    # parse the command line
    args = parser.parse_args()

    print args


    gz_pdb_files = get_pdb_files( pattern="ent.gz" )
    pdb_files = get_pdb_files( pattern=".pdb" )[0:]

    print "%s pdb files" % len( pdb_files )

    # print unzip_pdb( gz_pdb_files[0:10] )

    with Timer("get sstruc"):
        sstruc = sstruc_pdb( pdb_files )

    with Timer("struc records"):
        sstruc_records = [ (x[0], get_sstruc_records(x[1])) for x in sstruc ]

    with Timer("pickle records"):
        with open( "sstruc.pkl", "wb" ) as f:
            cPickle.dump( sstruc_records, f )

    # print sstruc_records

    with Timer("create sstruc db"):
        create_sstruc_db( sstruc_records )


    # print sstruc
    if False:
        for fpath, ss in sstruc_records:
            print fpath
            for elm in ss:
                print elm
    sstruc_count = sum([ len(x[1]) for x in sstruc ])
    print sstruc_count, sstruc_count/len(sstruc)


if __name__ == "__main__":
    main()