from __future__ import division

import collections

from utils.tool import _, _dir_init, CmdTool
from utils import numpdb
from utils.numpdb import NumPdb

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "dssp" )

DSSP_CMD = "mkdssp"



##################################################
# http://swift.cmbi.ru.nl/gv/dssp/index.html
##################################################
#
# COPYRIGHT
#   W. Kabsch, C. Sander and MPI-MF, 1983, 1985, 1988, 1994 1995
#   CMBI version by Elmar.Krieger@cmbi.kun.nl / November 18,2002
# USAGE
#   dsspcmbi [Options] PDB_File DSSP_File -> Read PDB_File and write DSSP_File
#   dsspcmbi [Options] -- dssp_file       -> Read from stdin and write DSSP_File
#   dsspcmbi -h                           -> Display this help screen
# OPTIONS
#   -na   Disables the calculation of accessible surface.
#   -c    Classic (old) format.
#   -w    Wide 2002 format (for future use,not the current standard).
#   -v    Verbose.
#   --    Read from standard input.

#   -h -? Prints a help message.
#   -V    Prints version, as in first line of the output.
# ADDITIONAL OPTIONS CONTRIBUTED BY DSSP USERS
#   By Emmanuel.Courcelle@toulouse.inra.fr
#   -ssa  Adds information about disulfide bonds to output file
#   -x    Renames residues with incomplete sidechains to 'X'
#   -alt2 Keeps an additional AltLoc indicator at the line ends


class Dssp( CmdTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    out = [
        _( "dssp_file", file="{pdb_file.stem}.dssp" )
    ]
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            DSSP_CMD, self.pdb_file, self.dssp_file
        ]
    def _post_exec( self ):
        npdb = NumPdb( self.pdb_file, features={
            "sstruc": False,
            "phi_psi": False,
            "extra": [( "dssp", "|S1", " " )]
        })
        self.get_numpdb_sstruc( npdb )
    def get_records( self ):
        return parse_dssp( self.dssp_file )
    def get_numpdb_sstruc( self, npdb ):
        return dssp2numpdb( self.get_records(), npdb )




# http://swift.cmbi.ru.nl/gv/dssp/HTML/descrip.html
#
# G: 3_10 helix
# H: alpha helix
# I: pi helix
# B: beta bridge
# E: beta sheet
# T: turn
# S: bend
#  : loop, coil
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
    return records


def dssp2numpdb( records, npdb ):
    d = collections.OrderedDict()
    for r in records:
        if r:
            d[( r[1], r[2] )] = r[3]
    for a in npdb._iter_resno():
        a['dssp'] = d.get(( a['resno'][0], a['chain'][0] ), ' ')
    return npdb


def dssp2pdb( records ):
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






