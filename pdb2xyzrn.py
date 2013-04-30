#! /usr/bin/env python


import re
import os

def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-calculate', help='do msms calculation', type=boolean)
    parser.add_argument(
        '-pdb_path', help='path to the pdb files and directories', type=str)
    parser.add_argument(
        '-test_sample', help='how many?', type=int)

    
    # parse the command line
    args = parser.parse_args()
    print args


    if args.calculate and args.pdb_path:

        pdb_files = get_pdb_files( args.pdb_path, pattern=".pdb" )
        if args.test_sample:
            pdb_files = pdb_files[0:args.test_sample]

        print "%s pdb files" % len( pdb_files )

        with Timer("calculate msms"):
            msms_ret = msms_pdb( pdb_files )



if __name__ == "__main__":
    main()