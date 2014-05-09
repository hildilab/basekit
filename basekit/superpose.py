from __future__ import with_statement
from __future__ import division



import os
import re
import json
import logging
import collections
import utils.numpdb as numpdb
import math
import shutil

import numpy as np
from numpy import *
from array import array

from basekit import utils
from utils import memoize_m, dir_walker, path
from utils.tool import _, _dir_init, CmdTool, ProviMixin, ParallelMixin
from utils.tool import RecordsMixin, PyTool
from utils.listing import merge_dic_list
from align import Muscle
import provi_prep as provi


pdbDIR, pdbPARENT_DIR, pdbTMPL_DIR = _dir_init( __file__, "pdb" )
DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "superpose" )
THESEUS_CMD = os.path.join( TMPL_DIR, "theseus" )

logging.basicConfig()
LOG = logging.getLogger('superpose')
LOG.setLevel( logging.ERROR )


class Theseus( CmdTool ):
    """A wrapper around the 'theseus' programm."""
    args = [
        _( "pdb_input", type="str", nargs="+", help="one directory (with at least 2 pdb-files) or at least 2 separate pdb-files" ),
        _( "dif_las|dls", type="bool", default=True, help="input files have different length and sequence" ),
        _( "prefix|pre", type="str", default="theseus" ),
        _( "amber", type="bool", default=False, help="if processing MD traces from AMBER" ),
        _( "selection|sa", type="int", range=[0,5], default=0, help="0(=alpha carbons for proteins, C1' atoms for nucleic acids - only option for superpositioning structures with different sequences), 1(=backbone), 2(=all), 3(=alpha and beta carbons), 4(=all heavy atoms - no hydrogens), default=0" ),
        _( "covariance|co", type="bool", default=False, help="use ML atomic covariance weighting (fit correlations, much slower), bad choice for many different structures with few residues fitting the correlation matrix, default=False" ),
        _( "embedding|e", type="int", default=2, help="embedding algorithm for initializing the average structure: 0(=none; use randomly chosen model), 2(=ML embedded structure), default=2" ),
        _( "first|f", type="bool", default=False, help="only read the first model of a multi-model PDB file, default=false" ),
        _( "hierarchical|g", type="int", default=1, help="hierarchical model for variances: 0(=none - may not converge), 1(=inverse gamma distribution), default=1" ),
        _( "iterations|i", type="int", default=200, help="maximum iterations, default=200" ),
        _( "constvar|k", type="int", default=-1, help="constant minimum variance, if set to negative value, the minimum varianve is determinded empirically, default=-1" ),
        _( "precision|p", type="str", default="1e-7", help="requested relative precision for convergence, default=1e-7" ),
        _( "ressele|s", type="str", default="all", help="Residue selection (e.g. A15-45:B50-55 = Residue 15 to 45 of chain A and residue 50-55 of chain B), default=all" ),
        _( "resexcl|S", type="str", default="none", help="Residue to exclude (e.g. A15-45:B50-55 = Residue 15 to 45 of chain A and residue 50-55 of chain B are excluded), default=none" ),
        _( "variance", type="bool", default=True, help="requested relative precision for convergence, default=True" ),
        _( "reference", type="file", default="", help="Reference file to superposition on, all rotations are relative to the first model in this file, default=first model" ),
        
    ]
    out = [
        _( "test_sup", file="test_sup.pdb" ),
        _( "test_sup_var", file="test_sup_var.pdb" ),
        _( "theseus_cmd", file="theseus_cmd.log"),
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        pdb_filess=[]
        all_pdb=[]
        if len(self.pdb_input)==1:
            v = ".*[pdb]$"
            for m, pdbfile in dir_walker( self.pdb_input[0], v ):
                stem=utils.path.stem(pdbfile)+utils.path.ext(pdbfile)
                #try:
                #    os.symlink(pdbfile, stem)
                #    shutil.move(stem, os.path.join(self.output_dir, stem))
                #except:
                #    pass
                pdb_filess.append(stem)
                all_pdb.append(pdbfile)
        else:
            for pdbfile in self.pdb_input:
                stem=utils.path.stem(pdbfile)+utils.path.ext(pdbfile)
                #try:
                #    os.symlink(pdbfile, stem)
                #    shutil.move(stem, os.path.join(self.output_dir, stem))
                #except:
                #    pass
                pdb_filess.append(stem)
                all_pdb.append(pdbfile)
        if self.dif_las:
            self.muscle = Muscle(
                self.prefix,
                pdb_files=all_pdb,
                mapfile=True,
                output_dir=self.output_dir
            )
            self.muscle()
            self.cmd = [ 
                THESEUS_CMD, 
                "-A", self.muscle.alignment_file, 
                "-M", self.muscle.map_file,
                "-a%d"%self.selection,
                "-e%d"%self.embedding,
                "-g%d"%self.hierarchical,
                "-i%d"%self.iterations,
                "-k%d"%self.constvar,
                "-p%s"%self.precision,
                "-r", self.prefix,
            ]
        else:
            self.cmd = [ 
                THESEUS_CMD, 
                "-a%d"%self.selection,
                "-e%d"%self.embedding,
                "-g%d"%self.hierarchical,
                "-i%d"%self.iterations,
                "-k%d"%self.constvar,
                "-p%s"%self.precision,
                "-r", self.prefix,
            ]
        if self.amber:
            self.cmd.append("--amber")
        if self.covariance:
            self.cmd.append("-c")
        if self.first:
            self.cmd.append("-f")
        if self.ressele!="all":
            self.cmd.append("-s%s"%self.ressele)
        if self.resexcl!="none":
            self.cmd.append("-S%s"%self.resexcl)
        if self.variance:
            self.cmd.append("-v")
        if self.reference!="":
            stem=self.reference.split("/")[-1]
            os.symlink(self.reference, stem)
            self.cmd=self.cmd+["-o", stem]
        v = ".*[pdb]$"
        self.cmd=self.cmd+pdb_filess
    def _post_exec( self ):
        if self.check( full=True ):
            newout=self.subdir("additional_files")
            for m, pdbfile in dir_walker( self.output_dir, "" ):
                stem=utils.path.stem(pdbfile)+utils.path.ext(pdbfile)
                if pdbfile.endswith("test_sup.pdb") or pdbfile.endswith("test_sup_var.pdb") or pdbfile.endswith("theseus_cmd.log"):
                    pass
                else:
                    shutil.move(pdbfile, os.path.join(newout, stem))



