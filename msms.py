#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division



import re
import os
import shutil
import argparse
import sqlite3
import functools
from string import Template

from utils import try_int, get_index, boolean
from utils.timer import Timer
from utils.job import run_command2, working_directory, do_parallel, _prep_func
from utils.db import get_pdb_files
from utils.jmol import run_jmol_script


MSMS_CMD = "msms"
BABEL_CMD = "babel"
TIMEOUT_CMD = "timeout"
PDB_TO_XYZR_CMD = "pdb_to_xyzr"




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


def pdb_to_xyzr( input_pdb, output_dir ):
    fdir, fname, pdb_id = pdb_path_split( input_pdb )
    xyzr = "%s.xyzr" % pdb_id
    pdb2 = "prep.pdb"
    cmd = [ TIMEOUT_CMD, "10", BABEL_CMD, '-i', 'pdb', pdb2, '-o', 'msms', xyzr ]
    with working_directory( output_dir ):
        pdb_select( input_pdb, pdb2 )
        run_command2( cmd, log="babel_cmd.log" )
    return xyzr


# def pdb_to_xyzr2( input_pdb, output_dir ):
#     fdir, fname, pdb_id = pdb_path_split( input_pdb )
#     out = "%s.xyzrn" % pdb_id

#     with open( "data/jmol/pdb2xyzrn.jspt", "r" ) as fp:
#         tpl_str = fp.read()
#     values_dict = {
#         "pdb_file": input_pdb,
#         "out_file": out
#     }
#     with working_directory( output_dir ):
#         run_jmol_script( Template( tpl_str ).substitute( **values_dict ) )
#     return out


def msms( input_pdb, output_dir ):
    xyzr = pdb_to_xyzr( input_pdb, output_dir )
    cmd = "%s -if %s -af area -of tri_surface" % (
        MSMS_CMD, xyzr
    )
    cmd2 = [ MSMS_CMD, "-if", xyzr, "-af", "area", "-of", "tri_surface" ]
    with working_directory( output_dir ):
        run_command2( cmd2, log="msms_cmd.log" )


@_prep_func
def _msms_pdb( fpath ):
    fdir, fname, pdb_id = pdb_path_split( fpath )
    output_dir = os.path.join( fdir, pdb_id, "msms" )
    if not os.path.exists( output_dir ):
        os.makedirs( output_dir )
    msms( os.path.abspath( fpath ), output_dir )
    return True

def msms_pdb( pdb_files, nworkers=None ):
    return do_parallel( _msms_pdb, pdb_files, nworkers )



def check_files( fpath ):
    fdir, fname, pdb_id = pdb_path_split( fpath )
    msms = os.path.exists( os.path.join( fdir, pdb_id, "msms" ) )
    xyzr = os.path.join( fdir, pdb_id, "msms", "%s.xyzr" % pdb_id )
    xyzr = os.path.exists( xyzr ) and os.path.getsize( xyzr )
    area = os.path.exists( os.path.join( fdir, pdb_id, "msms", "area.area" ) )
    return ( msms, xyzr, area )




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

    
    # parse the command line
    args = parser.parse_args()


    if args.pdb_path:
        pdb_files = get_pdb_files( args.pdb_path, pattern=".pdb" )
        if args.filter:
            pdb_filtered = []
            for fpath in pdb_files:
                msms, xyzr, area = check_files( fpath )
                if not msms or not area or not xyzr:
                    pdb_filtered.append( fpath )
            pdb_files = pdb_filtered
    if args.test_sample:
        pdb_files = pdb_files[0:args.test_sample]
        if args.test_sample<=10: print pdb_files
    if pdb_files:
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
            if args.status_out:
                with open( os.path.join( args.status_out, "calculate_msms.txt" ), "w" ) as fp:
                    fp.write( "\n".join([ "\t".join(map(str, x + check_files(x[0]))) for x in msms_ret ]) )


    if args.status and args.pdb_path:
        print "status"
        xyzr_fail = []
        msms_fail = []
        count_good = 0
        for fpath in pdb_files:
            msms, xyzr, area = check_files( fpath )
            if msms and area and xyzr:
                count_good += 1
            if msms and not xyzr:
                xyzr_fail.append( fpath )
            if msms and not area:
                msms_fail.append( fpath )

        print "calculated: %s, failed xyzr: %s, failed msms: %s" % ( 
            count_good, len( xyzr_fail ), len( msms_fail )
        )
        if args.status_out:
            with open( os.path.join( args.status_out, "failed_xyzr.txt" ), "w" ) as fp:
                fp.write( "\n".join( xyzr_fail ) )
            with open( os.path.join( args.status_out, "failed_msms.txt" ), "w" ) as fp:
                fp.write( "\n".join( msms_fail ) )


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