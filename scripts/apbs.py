#! /usr/bin/env python

import argparse
from basekit.apbs import Apbs





def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-pdb', help='pdb file', type=str)
    parser.add_argument(
        '-output_dir', help='output dir', type=str, default=".")


    # parse the command line
    args = parser.parse_args()

    if args.pdb and args.output_dir:
    	Apbs( args.pdb, output_dir=args.output_dir )



if __name__ == "__main__":
    main()