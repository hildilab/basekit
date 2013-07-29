from __future__ import division

import utils.path
from utils.tool import _, _dir_init, CmdTool

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "msa" )

MUSCLE_CMD = "muscle"


# MUSCLE v3.8.31 by Robert C. Edgar

# http://www.drive5.com/muscle
# This software is donated to the public domain.
# Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.


# Basic usage

#     muscle -in <inputfile> -out <outputfile>

# Common options (for a complete list please see the User Guide):

#     -in <inputfile>    Input file in FASTA format (default stdin)
#     -out <outputfile>  Output alignment in FASTA format (default stdout)
#     -diags             Find diagonals (faster for similar sequences)
#     -maxiters <n>      Maximum number of iterations (integer, default 16)
#     -maxhours <h>      Maximum time to iterate in hours (default no limit)
#     -html              Write output in HTML format (default FASTA)
#     -msf               Write output in GCG MSF format (default FASTA)
#     -clw               Write output in CLUSTALW format (default FASTA)
#     -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
#     -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
#     -quiet             Do not write progress messages to stderr
#     -version           Display version information and exit



class Muscle( CmdTool ):
    args = [
        _( "fasta_file", type="file", ext="fasta" )
    ]
    out = [
        _( "muscle_file", file="muscle_{fasta_file.stem}.fasta" )
    ]
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            MUSCLE_CMD, "-in", self.fasta_file, "-out", self.muscle_file
        ]

