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
from utils.bio import AA1
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
        _( "chain|ch", type="str", default="all", help="for each one specific (e.g.: AB* for 3 files, first one chain A, second chain B and third one all)" ),
        _( "amber", type="bool", default=False, help="if processing MD traces from AMBER" ),
        _( "selection|sa", type="int", range=[0,5], default=0, help="0(=alpha carbons for proteins, C1' atoms for nucleic acids - only option for superpositioning structures with different sequences), 1(=backbone), 2(=all), 3(=alpha and beta carbons), 4(=all heavy atoms - no hydrogens), default=0" ),
        _( "covariance|co", type="bool", default=False, help="use ML atomic covariance weighting (fit correlations, much slower), bad choice for many different structures with few residues fitting the correlation matrix, default=False" ),
        _( "embedding|e", type="int", default=2, help="embedding algorithm for initializing the average structure: 0(=none; use randomly chosen model), 2(=ML embedded structure), default=2" ),
        _( "first|f", type="bool", default=True, help="option only for input with same length and sequence: only read the first model of a multi-model PDB file, default=True" ),
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
        _( "test_sup", file="{prefix}_sup.pdb" ),
        _( "test_sup_var", file="{prefix}_sup_var.pdb" ),
        _( "theseus_cmd", file="theseus_cmd.log"),
        _( "variances_file", file="{prefix}_variances.txt" ),
        _( "out_test_sup", file="{prefix}_var.pdb" ),
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        self.stem_ini()
        if self.dif_las:
            self.muscle = Muscle(
                self.prefix,
                pdb_files=self.pdb_input,
                mapfile=True,
                chain=self.chain,
                output_dir=self.outpath( "muscle" ),
                run=False
            )
            self.output_files += self.muscle.output_files
            self.sub_tool_list.append( self.muscle )
            
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
                "-O", "-P 6",
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
                "-O", "-P 6",
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
            self.sstem=self.reference.split("/")[-1]
            os.symlink(self.reference, self.sstem)
        else:
            self.sstem=self.stem_list[0]
        self.cmd=self.cmd+["-o", self.sstem]
        self.cmd=self.cmd+self.stem_list
    def _pre_exec( self ):
        if self.dif_las:
            self.muscle()
            self.pdb_input=self.muscle.new_files
        for pdbfile in self.pdb_input:
            name=utils.path.stem(pdbfile)+utils.path.ext(pdbfile)
            try:
                os.symlink(pdbfile, name)
                shutil.move(name, os.path.join(self.output_dir, name))
            except:
                pass
    def _post_exec( self ):
        #generates a file with the variances in the bfactor-field
        if self.dif_las:
            #writes variances into a vector
            variances=[]
            with open(self.variances_file, "r") as fp:
                for line in fp:
                    if line.startswith("RES"):
                        line=line.split()
                        variances.append([line[2],line[3],line[4]])
            #gets the alignments needed for matching the variances to the residues
            aligns=[]
            ids=[]
            with open(self.muscle.alignment_file, "r") as fp:
                counter=0
                alin=[]
                for line in fp:
                    if line.startswith('>'):
                        fi=str.strip(line).split('>')[1]
                        if fi==self.sstem:
                            pos_leading=counter
                        ids.append(fi)
                        if counter!=0:
                            aligns.append(alin)
                        counter=counter+1
                        alin=[]
                    else:
                        for elem in str.strip(line):
                            alin.append(elem)
                aligns.append(alin)
            #generates a list containing a list and id for each file:
            # the list contains the var for each residue, depending,
            # if in the alignment of the leading structure (self.sstrem)
            # '-' or a residue is given, because there are only variances
            # for the leading structure
            
            #open question: how to set '-'? for very unequal files the variances
            # are very small - problems with comparison?
            rightordervarlist=[]
            for fi in range(len(ids)):
                position=0
                i=0
                rightordervar=[]
                for index, residue in enumerate(aligns[pos_leading]):
                    if residue=='-':
                        if aligns[fi][index]!='-':
                            rightordervar.append('    -1')
                    else:
                        if variances[position][1]==1:
                            i=1
                        if aligns[fi][index]!='-':
                            if str(i)==variances[position][1]:
                                var="%.2f" % float(variances[position][2])
                                if len(var)<6:
                                    leerz=" "*(6-len(var))
                                    var=leerz+var
                                rightordervar.append(var)
                                position+=1
                            else:
                                rightordervar.append('    -1')
                        i+=1
                rightordervarlist.append([rightordervar, ids[fi]])
            #writes the variances into the *var.pdb-file
            with open(self.test_sup, "r") as fp, open(self.out_test_sup, "w") as fp2:
                nums=[]
                for lin in fp:
                    if lin.startswith('ATOM'):
                        act=lin[22:26]
                        if not(old=='' or act==old):
                            count+=1
                        old=act
                        var=mod_variances[count]
                        line=lin[:60]+var+lin[66:]
                        fp2.write(line)
                    elif lin.startswith('MODEL'):
                        model=nums[int(lin.split()[1])-1]
                        for mod in rightordervarlist:
                            if model==mod[1]:
                                mod_variances=mod[0]
                        count=0
                        old=''
                        fp2.write(lin)
                    elif lin.startswith('REMARK   MODEL'):
                        nums.append(lin.split()[3])
                        fp2.write(lin)
                    else:
                        fp2.write(lin)
            #sorts files into additional dir
            newout=self.subdir("additional_files")
            for m, pdbfile in dir_walker( self.output_dir, "" ):
                stem=utils.path.stem(pdbfile)+utils.path.ext(pdbfile)
                if "muscle" in pdbfile or pdbfile.endswith(self.prefix+"_var.pdb") or pdbfile.endswith(self.prefix+"_sup.pdb") or pdbfile.endswith(self.prefix+"_sup_var.pdb") or pdbfile.endswith("theseus_cmd.log") or pdbfile.endswith("_variances.txt"):
                    pass
                else:
                    shutil.move(pdbfile, os.path.join(newout, stem))
    def stem_ini(self):
        #gets all stems of all files
        self.stem_list=[]
        if len(self.pdb_input)==1:
            v = ".*[pdb]$"
            for m, pdbfile in dir_walker( self.pdb_input[0], v ):
                name=utils.path.stem(pdbfile)+utils.path.ext(pdbfile)
                self.stem_list.append(name)
        elif isinstance(self.pdb_input, list):
            for pdbfile in self.pdb_input:
                name=utils.path.stem(pdbfile)+utils.path.ext(pdbfile)
                self.stem_list.append(name)
    



