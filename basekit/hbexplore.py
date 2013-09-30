from __future__ import with_statement
from __future__ import division


import collections


from utils import memoize_m
from utils.tool import _, _dir_init, CmdTool, ProviMixin

import provi_prep as provi


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "hbexplore" )
HBEXPLORE_CMD = "HBX"



# 1=donor, 2=acceptor ???
HbxHbond = collections.namedtuple( "HbxHbond", [
    "chain1", "resno1", "resname1", "atomname1",
    "chain2", "resno2", "resname2", "atomname2",
    "dDA", "dHA", "aDHA", "aHAhyb", "type"
])

def parse_hbx_output( hbexplore_file ):
    """Function that parses HBexplores .anal files"""
    hbonds = []
    id_dict = {}
    hbpart=0
    #   SG  CYS   127      O   ILE   124    3.2597     
    # TOM     18  CB  VAL     2       3.462  14.126   6.256  1.00 14.07      193L 162
    # parse HBX.anal file
    #87     NH2 ARG A  71      O   MET A   1    3.0821  2.356   128.0    53.8   s-s
    with open( hbexplore_file, "r" ) as fp:
        for l in fp:
            if not l:
                continue
            if l.startswith('No.'):
                hbpart=1
            elif hbpart and l[0]=='=':
                break
            elif not hbpart or len(l)<=30:
                continue
            elif hbpart and l[0]!='-':
                # print l[:-1]
                ts = l.split("\t")
                t=ts[1]
                hb_id = int(ts[0])
                if( hb_id not in id_dict ):
                    id_dict[ hb_id ] = True
                    hbonds.append( HbxHbond(
                        t[12:14].strip(),   # chain1
                        int(t[14:18]),      # resno1
                        t[8:11].strip(),    # resname1
                        t[4:7].strip(),     # atomname1
                        t[31:33].strip(),   # chain2
                        int(t[33:37]),      # resno2
                        t[27:31].strip(),   # resname2
                        t[23:27].strip(),   # atomname2
                        float(t[40:]),      # dDA
                        float(ts[2]),       # dHA
                        float(ts[3]),       # aDHA
                        float(ts[4]),       # aHAhyb
                        ts[5][:-1]          # type
                    ))
    return hbonds



class HBexplore( CmdTool, ProviMixin ):
    """A wrapper around the 'hbexplore' programm."""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "DAmax", type="float", default=3.9 ),
        _( "HAmax", type="float", default=2.5 ),
        _( "DHAmin", type="int", default=90 ),
        _( "HAA1_DAA1min", type="int", default=90 ),
        _( "HAhybmax", type="int", default=60 ),
        _( "waterAsDonorAcceptor", type="int", default=1 ),
        _( "c_hAsDonor", type="int", default=2 ),
    ]
    out = [
        _( "params_file", file="{pdb_file.stem}_hbx.params" ),
        _( "hbx_file", file="{pdb_file.basename}_hbx.anal" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "hbexplore.provi"
    def _init( self, *args, **kwargs ):
        self.cmd = [ HBEXPLORE_CMD ]
        self.cmd_input = self.params_file
        self.no_cmd = False
    def _pre_exec( self ):
        with open(self.params_file, "w") as fp:
            fp.write( "\n".join([
                self.relpath( self.pdb_file ),
                ".",            # rel output dir
                "n",            # use_default_params
                str(self.DAmax),
                str(self.HAmax),
                str(self.DHAmin),
                str(self.HAA1_DAA1min),
                str(self.HAhybmax),
                str(self.waterAsDonorAcceptor),
                str(self.c_hAsDonor),
                "",             # separate_tables_and_distributions,
                "0",            # additional_connectivity_files,
                "",             # potential_hydrogen_bond_dataset,
                "n",            # search_for_specific_hbonds
            ]) + "\n\n" )
    def _post_exec( self ):
        provi.prep_hbexplore( self.hbx_file, self.pdb_file )
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            hbonds_file=self.relpath( self.hbx_file + ".bonds" )
        )
    @memoize_m
    def get_hbonds( self ):
        return parse_hbx_output( self.hbx_file )



