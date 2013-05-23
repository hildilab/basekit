from __future__ import with_statement
from __future__ import division



import re
import os
import shutil
import argparse

from utils.tool import CmdTool


DSSP_CMD = "dssp"



##################################################
# http://swift.cmbi.ru.nl/gv/dssp/index.html
##################################################
#
#   dssp [-na] [-v] pdb_file [dssp_file]
#   dssp [-na] [-v] -- [dssp_file]
#   dssp [-h] [-?] [-V]
#
# The possible DSSPold commandline options are:

# -na     Disables the calculation of accessible surface.
# -c      Classic (pre-July 1995) format.
# -v      Verbose.
# --      Read from standard input.
# -h      Prints a help message.
# -?      Same as -h
# -l      Prints the license information.
# -V      Prints version, as in first line of the output.


class Dssp( CmdTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" }
    ]
    def _init( self, pdb_file, **kwargs ):
    	self.pdb_file = os.path.abspath( pdb_file )
    	stem = os.path.splitext( os.path.split( self.pdb_file )[-1] )[0]
        self.dssp_file = "%s.dssp" % stem
        self.cmd = [ 
            DSSP_CMD, self.pdb_file, self.dssp_file
        ]
        self.output_files = [ self.dssp_file ]






# http://swift.cmbi.ru.nl/gv/dssp/HTML/descrip.html
#
def parse_dssp( dssp_path ):
    records = [ None ]
    with open( dssp_path, "r" ) as fp:
        residue_section = False
        for line in fp:
            if line.startswith( "  #  RESIDUE AA STRUCTURE" ):
                residue_section = True
            elif line.strip() and residue_section:
                if line[13]=="!":
                    records.append( None )
                else:
                    records.append([
                        int(line[0:5]),                 # no
                        int(line[5:10]),                # resno
                        line[11],                       # chain
                        line[16],                       # ss type
                        int(line[25:29]),               # hbond to resno
                        line[13],                       # aa1
                        line[33],                       # beta sheet label
                    ])
    prev_ss = " "
    prev_chain = " "
    sstruc = [[]]
    for i, r in enumerate( records ):
        if r:
            if r[3]==prev_ss:
                sstruc[-1].append(r)
            else:
                sstruc.append([r])
            prev_chain = r[2]
            prev_ss = r[3]
    sheet_label_count = defaultdict(int)
    for ss in sstruc:
        if ss[0][3] in [ "E" ]:
            sheet_label_count[ ss[0][6] ] += 1
    helix_counter = 0
    helix_list = []
    sheet_counter = 0
    sheet_list = []
    for ss in sstruc:
        if ss[0][3] in [ "H" ]:
            helix_counter += 1
            helix_list.append(
                "HELIX% 5i% 4i %s %s% 5i  %s %s% 5i % 2i                               % 5i    " % (
                    helix_counter, helix_counter, 
                    numpdb.AA3.get( ss[0][5], "XXX" ), ss[0][2], ss[0][1], 
                    numpdb.AA3.get( ss[-1][5], "XXX" ), ss[-1][2], ss[-1][1],
                    1, len(ss)
                )
            )
        if ss[0][3] in [ "E" ]:
            sheet_counter += 1
            e = records[ ss[0][4] ]
            e2 = records[ ss[-1][4] ]
            if not e2 or ss[0][0]<e2[0]:
                hbond = "                                        "
            else:
                hbond = "  O  %s %s% 4i   N  %s %s% 4i           " % (
                    numpdb.AA3.get( ss[-1][5], "XXX" ), ss[-1][2], ss[-1][1],
                    numpdb.AA3.get( e2[5], "XXX" ), e2[2], e2[1],
                )
            sheet_list.append(
                "SHEET% 5i   %s% 2i %s %s% 4i  %s %s% 4i -1%s" % (
                    sheet_counter, ss[0][6], sheet_label_count[ ss[0][6] ],
                    numpdb.AA3.get( ss[0][5], "XXX" ), ss[0][2], ss[0][1], 
                    numpdb.AA3.get( ss[-1][5], "XXX" ), ss[-1][2], ss[-1][1],
                    hbond
                )
            )
    # print "HELIX    1   1 GLU A   33  HIS A   65  1                                  33    "
    # for h in helix_list:
    #     print h
    # print "SHEET    2   A 2 TYR A  10  VAL A  11 -1  O  VAL A  11   N  THR A   4           "
    # for e in sheet_list:
    #     print e
    return "\n".join( helix_list + sheet_list )






