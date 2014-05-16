from __future__ import with_statement
from __future__ import division

import os
import logging
import collections

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
    out = [
        _( "theseus_cmd", file="theseus_cmd.log"),
    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            THESEUS_CMD, 
            "-f", "-F",
            ]
        self.cmd=self.cmd+self.pdb_input


def pdblist2fasta(inputfiles, fasta_file, outdir, mapfile=""):
    if mapfile!="":
        fp2=open(mapfile, "w")
    stem_list=[]
    for index, fi in enumerate(inputfiles):
        stem=utils.path.stem(fi)
        name=stem+utils.path.ext(fi)
        if mapfile!="":
            fp2.write(fi.split("/")[-1]+" "+name+"\n")#stem+"\n")
        npdb = numpdb.NumPdb( fi )
        new_pdb=os.path.join(outdir, name)
        new=npdb.copy(**{"record":"ATOM  "})
        new.write(new_pdb)
        sequence=new.sequence()
        sequence=sequence.replace("?","")
        identifier=">%s\n"%name#stem
        stem_list.append(name)
    if mapfile!="":
        fp2.close()
    fasta = TheseusMakeFasta(
        stem_list,
        output_dir=outdir
    )
    fasta()
    with open(fasta_file, "w") as fp:
        v = ".*[fst]$"
        for m, pdbfile in dir_walker( outdir, v ):
            with open(pdbfile, "r") as fp3:
                content=fp3.read()
                fp.write(content)
                fp.write("\n")


class Muscle( CmdTool):
    """A wrapper around the 'Muscle' programm."""
    args = [
        _( "fasta_file", type="file", ext="fasta", help="FASTA-file with all sequences of the to be aligned pdb-files OR if using pdb_files: name of the file" ),
        _( "pdb_files|pf", type="str", nargs="+", default="", help="at least 2 pdb-files or a directory containing at least 2 pdb-files to be aligned" ),
        _( "diags", type="bool", default=False, help="Find diagonals (faster for similar sequences), (bool, default=False)" ),
        _( "maxiters", type="int", default=16, help="Maximum number of iterations (integer, default=16)" ),
        _( "maxhours", type="int", default=0, help="Maximum time to iterate in hours (default no limit)" ),
        _( "format", type="str", options=["fasta", "html", "msf", "clw"], default="fasta", help="[fasta, html, msf(GCG MSF), clw(CLUSTALW)] (default fasta)" ),
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
        if self.pdb_files!="":
            files=[]
            if len(self.pdb_files)==1:
                v = ".*[pdb]$"
                for m, pdbfile in dir_walker( self.pdb_files, v ):
                    files.append(pdbfile)
            else:
                files=self.pdb_files
            if self.mapfile:
                pdblist2fasta(files, self.fasta_file_out, self.output_dir, self.map_file)
            else:
                pdblist2fasta(files, self.fasta_file_out, self.output_dir)
            self.fasta_file=self.fasta_file_out
        if self.diags:
            if self.maxhours!=0:
                self.cmd = [ 
                    MUSCLE_CMD,
                    "-in", self.fasta_file,
                    "-out", self.alignment_file,
                    "-diags"
                    "-maxiters", self.maxiters,
                    "-maxhours", self.maxhours,
                    "-%s"%self.format,
                    "-log", self.log_file2
                ]
            else:
                self.cmd = [ 
                    MUSCLE_CMD,
                    "-in", self.fasta_file,
                    "-out", self.alignment_file,
                    "-diags"
                    "-maxiters", self.maxiters,
                    "-%s"%self.format,
                    "-log", self.log_file2
                ]
        else:
            if self.maxhours!=0:
                self.cmd = [ 
                    MUSCLE_CMD,
                    "-in", self.fasta_file,
                    "-out", self.alignment_file,
                    "-maxiters", self.maxiters,
                    "-maxhours", self.maxhours,
                    "-%s"%self.format,
                    "-log", self.log_file2
                ]
            else:
                self.cmd = [ 
                    MUSCLE_CMD,
                    "-in", self.fasta_file,
                    "-out", self.alignment_file,
                    "-maxiters", self.maxiters,
                    "-%s"%self.format,
                    "-log", self.log_file2
                ]



