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
from utils.job import _prep_func, do_parallel
from utils.math import vec_angle

import utils.numpdb as numpdb






def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-pdb', help='analyze a pdb file', type=str)

    
    # parse the command line
    args = parser.parse_args()



    if args.pdb:
        npdb = numpdb.NumPdb( args.pdb )
        print npdb.dist( { "chain": "A" }, { "chain": "B" } )


if __name__ == "__main__":
    main()