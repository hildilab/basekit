#! /usr/bin/env python

import argparse
from basekit.spider import Spider, SpiderConvert, SpiderDeleteFilledDensities, SpiderBox, SpiderCrosscorrelation





def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-mrc', help='mrc file', type=str)
    parser.add_argument(
        '-map', help='map file', type=str)
    parser.add_argument(
        '-pdb', help='pdb file', type=str)
    parser.add_argument(
        '-res1', help='res1', type=str)
    parser.add_argument(
        '-res2', help='res2', type=str)
    parser.add_argument(
        '-len', help='loop length', type=int)
    parser.add_argument(
        '-map1', help='map1 file', type=str)
    parser.add_argument(
        '-map2', help='map2 file', type=str)
    parser.add_argument(
        '-box', help='box variables file', type=str)
    parser.add_argument(
        '-loops', help='loops directory', type=str)
    parser.add_argument(
        '-o', '--output_dir', help='output dir', type=str, default=".")


    # parse the command line
    args = parser.parse_args()

    if args.mrc and args.output_dir:
    	SpiderConvert( args.mrc, output_dir=args.output_dir )

    if args.map and args.pdb and not args.res1 and args.output_dir:
        SpiderDeleteFilledDensities( 
            args.map, args.pdb, output_dir=args.output_dir
        )

    if args.map and args.pdb and args.res1 and args.res2 and args.len and args.output_dir:
        SpiderBox( 
            args.map, args.pdb, args.res1, args.res2, args.len,
            output_dir=args.output_dir
        )

    if args.map1 and args.map2 and args.box and args.loops:
        SpiderCrosscorrelation(
            args.map1, args.map2, args.box, args.loops,
            output_dir=args.output_dir
        )


if __name__ == "__main__":
    main()