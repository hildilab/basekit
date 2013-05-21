#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division


import re
import sys
import os
import shutil
import argparse
import sqlite3
import functools
from itertools import izip
from operator import itemgetter, methodcaller
from string import Template

from utils import try_int, get_index, boolean, working_directory, copy_dict
from utils.timer import Timer
from utils.job import do_parallel
from utils.db import get_pdb_files
from utils.tool import CmdTool, make_args


MSMS_CMD = "msms"
BABEL_CMD = "babel"



def pdb_select( input_pdb, output_pdb ):
    with open( input_pdb, "r" ) as fp_in:
        with open( output_pdb, "w" ) as fp_out:
            for line in fp_in:
                if line.startswith("ATOM"):
                    fp_out.write( line )
                elif line.startswith("ENDMDL"):
                    break


def pdb_path_split( fpdb ):
    fdir, fname = os.path.split( fpdb )
    pdb_id = fname.split(".")[0]
    return ( fdir, fname, pdb_id )


class Msms( CmdTool ):
    args = make_args([
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "density", "type": "slider", "range": [1, 10], "fixed": True, "default_value": 1.0  }
    ])
    def _init( self, pdb_file, density=1.0, **kwargs ):
        self.pdb2xyzr = Pdb2xyzr( pdb_file, **copy_dict( kwargs, run=False ) )
        self.cmd = [ 
            MSMS_CMD, "-if", self.pdb2xyzr.xyzr_file, 
            "-af", "area", "-of", "tri_surface", "-density", str(density)
        ]
        self.output_files = self.pdb2xyzr.output_files + \
            [ "area.area", "tri_surface.face", "tri_surface.vert" ]
    def _pre_exec( self ):
        self.pdb2xyzr()



class Pdb2xyzr( CmdTool ):
    args = make_args([
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ])
    def _init( self, pdb_file, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        self.pdb_prep_file = "prep.pdb"
        self.xyzr_file = "%s.xyzr" % os.path.splitext( os.path.split( self.pdb_file )[-1] )[0]
        self.cmd = [ 
            BABEL_CMD, '-i', 'pdb', self.pdb_prep_file,
            '-o', 'msms', self.xyzr_file 
        ]
        self.output_files = [ self.pdb_prep_file, self.xyzr_file ]
    def _pre_exec( self ):
        with working_directory( self.output_dir ):
            pdb_select( self.pdb_file, self.pdb_prep_file )



def _msms_pdb( fpath, run=True ):
    fdir, fname, pdb_id = pdb_path_split( fpath )
    output_dir = os.path.join( fdir, pdb_id, "msms" )
    return Msms( fpath, output_dir=output_dir, run=run, timeout=10 )


def msms_pdb( pdb_files, nworkers=None, run=True ):
    return do_parallel( _msms_pdb, pdb_files, nworkers, run=run )




def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-calculate', help='do msms calculation', type=boolean)
    parser.add_argument(
        '-status', help='summary of how many files are already calculated', type=boolean)
    parser.add_argument(
        '-status_out', help='output dir for status files', type=str)
    parser.add_argument(
        '-analyze', help='do msms analysis', type=boolean)
    parser.add_argument(
        '-filter', help='filter calculated', type=boolean)
    parser.add_argument(
        '-delete', help='delete msms directory', type=boolean)
    parser.add_argument(
        '-pdb_path', help='path to the pdb files and directories', type=str)
    parser.add_argument(
        '-test_sample', help='how many?', type=int)
    parser.add_argument(
        '-pdb', help='analyze a single pdb file', type=str)
    parser.add_argument(
        '-out', help='output directory', type=str, default=".")

    
    # parse the command line
    args = parser.parse_args()

    if args.pdb and args.out:
        with Timer("single pdb"):
            # Msms( args.pdb, output_dir=args.out )
            Msms( fileargs=True, output_dir=args.out )

    if args.pdb_path:
        pdb_files = get_pdb_files( args.pdb_path, pattern=".pdb" )
        if args.filter:
            pdb_filtered = []
            for fpath, d in zip( pdb_files, msms_pdb( pdb_files, run=False ) ):
                if not d.check(): pdb_filtered.append( fpath )
            pdb_files = pdb_filtered
    if args.test_sample:
        pdb_files = pdb_files[0:args.test_sample]
        if args.test_sample<=10: print pdb_files
    if "pdb_files" in locals():
        print "%s pdb files" % len( pdb_files )


    if args.delete and args.pdb_path:
        for fpath in pdb_files:
            fdir, fname, pdb_id = pdb_path_split( fpath )
            msms_dir = os.path.join( fdir, pdb_id, "msms" )
            if os.path.exists( msms_dir ):
                shutil.rmtree( msms_dir )


    if args.calculate and args.pdb_path:
        with Timer("calculate msms"):
            msms_ret = msms_pdb( pdb_files )


    if args.status and args.pdb_path:
        print "[ status ]"

        data_list = msms_pdb( pdb_files, run=False )        
        of = data_list[0].output_files
        check_counter = [0]*( len(of)+1 )
        if args.status_out:
            check_list = []
        for d in data_list:
            check = d.check( full=True )
            for i, f in enumerate(of):
                if not check[i]: check_counter[i+1] += 1
            if not all( check ):
                check_counter[0] += 1
                if args.status_out:
                    check_list.append( 
                        "%s\t%s" % ( 
                            d.output_dir, 
                            "\t".join( map(str, map(int, check)) ) 
                        )
                    )

        print "  failed: %i" % check_counter[0]
        for i, f in enumerate(of):
            print "  failed: [%s]: %i" % ( f, check_counter[i+1] )

        if args.status_out:
            with open( os.path.join( args.status_out, "msms_failed.txt" ), "w" ) as fp:
                fp.write( "path\t%s\n" % "\t".join( of ) )
                fp.write( "\n".join( check_list ) )



    if args.analyze and args.pdb_path:
        for fpath in pdb_files:
            fdir, fname, pdb_id = pdb_path_split( fpath )
            area_file = os.path.join( fdir, pdb_id, "msms", "area.area" )
            if os.path.exists( area_file ):
                with open( area_file ) as fp:
                    for line in fp:

                        print line



if __name__ == "__main__":
    main()