from __future__ import with_statement
from __future__ import division

import os
import logging
import collections
import shutil

from basekit import utils
import utils.numpdb as numpdb
from utils import dir_walker
from utils.tool import _, _dir_init, CmdTool

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "superpose" )
THESEUS_CMD = os.path.join( TMPL_DIR, "theseus" )

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "align" )
MUSCLE_CMD = os.path.join( TMPL_DIR, "muscle3.8.31" )


logging.basicConfig()
LOG = logging.getLogger('align')
LOG.setLevel( logging.ERROR )



class TheseusMakeFasta( CmdTool ):
    """A wrapper around the 'theseus-align' programm."""
    args = [
        _( "pdb_input", type="str", nargs="+", help="list of pdb-files" ),
    ]
    out = []
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        self.pdb_input2 = []
        for fi in self.pdb_input:
            name=utils.path.stem(fi)+utils.path.ext(fi)
            try:
                os.symlink(fi, name)
                shutil.move(name, os.path.join(self.output_dir, name))
            except:
                pass
            self.pdb_input2.append(name)
        self.cmd = [ 
            THESEUS_CMD, 
            "-f", "-F",
            ]
        self.cmd=self.cmd+self.pdb_input2


class Muscle( CmdTool):
    """A wrapper around the 'Muscle' programm."""
    args = [
        _( "fasta_file", type="file", help="FASTA-file with all sequences of the to be aligned pdb-files OR if using pdb_files: name of the file" ),
        _( "pdb_files|pf", type="str", nargs="+", default="", help="at least 2 pdb-files or a directory containing at least 2 pdb-files to be aligned" ),
        _( "diags", type="bool", default=False, help="Find diagonals (faster for similar sequences), (bool, default=False)" ),
        _( "maxiters", type="int", default=16, help="Maximum number of iterations (integer, default=16)" ),
        _( "maxhours", type="int", default=0, help="Maximum time to iterate in hours (default no limit)" ),
        _( "format", type="str", options=["fasta", "html", "msf", "clw"], default="fasta", help="[fasta, html, msf(GCG MSF), clw(CLUSTALW)] (default fasta)" ),
        _( "chain|ch", type="str", default="all", help="for each one specific (e.g.: AB* for 3 files, first one chain A, second chain B and third one all)" ),
        _( "mapfile", type="bool", default=False, help="generates a 'map' file for mapping inside Theseus(default=False)" ),
    ]
    out = [
        _( "alignment_file", file="{fasta_file.stem}_muscle.aln" ),
        _( "fasta_file_out", file="{fasta_file.stem}.fasta", optional=True  ),
        _( "log_file", file="muscle_cmd.log" ),
        _( "log_file2", file="{fasta_file.stem}_muscle.log" ),
        _( "map_file", file="{fasta_file.stem}_muscle.filemap", optional=True ),

    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        #generates the Fasta-File with Theseus
        if self.pdb_files!="":
            self.prep()
            #sets up TheseusMakeFasta
            self.fasta = TheseusMakeFasta(
                self.new_files,
                output_dir=self.outpath( "fasta" ),
                run=False
            )
            self.output_files += self.fasta.output_files
            self.fasta_file=self.fasta_file_out
            self.sub_tool_list.append( self.fasta )

        self.cmd = [ 
            MUSCLE_CMD,
            "-in", self.fasta_file,
            "-out", self.alignment_file,
            "-maxiters", self.maxiters,
            "-%s"%self.format,
            "-log", self.log_file2
        ]
        if self.diags:
            self.cmd.append("-diags")
        if self.maxhours!=0:
            self.cmd.append("-maxhours")
            self.cmd.append(self.maxhours)
        
    def prep( self ):
        files=[]
        #get all files from dir
        if len(self.pdb_files)==1:
            v = ".*[pdb]$"
            for m, pdbfile in dir_walker( self.pdb_files, v ):
                files.append(pdbfile)
        else:
            files=self.pdb_files
        if self.mapfile:
            fp2=open(self.map_file, "w")
        self.new_files=[]
        #generating a new pdb with only ATOM and if wanted selected chains
        for index, fi in enumerate(self.pdb_files):
            name=utils.path.stem(fi)+utils.path.ext(fi)
            if self.mapfile:
                fp2.write(fi.split("/")[-1]+" "+name+"\n")
            npdb = numpdb.NumPdb( fi )
            new_pdb=os.path.join(self.output_dir, name)
            if self.chain!="all" and self.chain[index]!="*":
                new=npdb.copy(**{"record":"ATOM  ", "chain":self.chain[index]})
            else:
                new=npdb.copy(**{"record":"ATOM  "})
            new.write(new_pdb)
            self.new_files.append(new_pdb)
        if self.mapfile:
            fp2.close()
    def _pre_exec( self ):
        if self.pdb_files!="":
            self.fasta()
            #merges all fasta files
            with open(self.fasta_file_out, "w") as fp:
                v = ".*[fst]$"
                for m, pdbfile in dir_walker( self.output_dir, v ):
                    with open(pdbfile, "r") as fp3:
                        content=fp3.read()
                        fp.write(content)
                        fp.write("\n")


