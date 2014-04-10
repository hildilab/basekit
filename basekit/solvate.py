import os

from utils.tool import _, _dir_init, CmdTool


pdbDIR, pdbPARENT_DIR, pdbTMPL_DIR = _dir_init( __file__, "pdb" )
DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "solvate" )
solvate_CMD = os.path.join( TMPL_DIR, "solvate" )


class Solvate( CmdTool):
    """
    SOLVATE is a program to construct an atomic solvent environment model
    for a given atomic macromolecule model (solute) for use in molecular
    dynamics simulations.
        """
    args = [
        _( "pdb_input", type="file",
            help="specifies the name of the input pdb-/psf-file of the solute. "+
                 "If the -ion option is set, both pdb- and psf-file are required; "+
                 "otherwise only the pdb-file is needed. No extension must be given! "+
                 "The input file name may be omitted, in which case SOLVATE creates "+
                 "a pure spherical water droplet centered at the origin." ),
        _( "thick|th", type="float", default=10.0,
            help="specifies the minimum water shell thickness in Angstrom (dafult: 10.0),"+
                 "Nowhere will the solute be closer to the surface of the solvent than this value"+
                 "Attention: Values smaller than 3.0 may confuse the program." ),
        _( "radius|r", type="float", default=100000.0,
            help="specifies the maximum boundary curvature radius of the solvent surface "+
                 "in Angstrom (dafault: 100000.0). A very large value (e.g. default) means "+
                 "that the surface can have rather flat parts. Smaller values generate 'rounder' "+
                 "solvent volumes, and, as a consequence, more solvent molecules. The value "+
                 "should not be considerably smaller than the size of the solute." ),
        _( "ngauss|n", type="int", default=1,
            help="specifies the number of gaussians to be used to define the solvent "+
                 "surface (default: True). The more detailed irregular features of the "+
                 "surface are required, the more gaussians should be used. Consequently, "+
                 "thin water shells typically require more gaussians than thicker ones. "+
                 "Note, however, that a large number may slow down subsequent molecular "+
                 "dynamics simulations. Typical values for ngauss are five to ten" ),
        _( "ug", type="bool", default=False,
            help="If the -ug (use gaussians) option is set, all steps required for "+
                 "the computation of the approximate density function f (STEPs 2 and 3) "+
                 "are skipped, and the required parameters for defining f are read from "+
                 "the file gaussians.lis instead, which is always written when f is computed. "+
                 "The reason for this option is that STEP 3 is quite time consuming." ),
        _( "ub", type="bool", default=False,
            help="If the -ub (use boundary) option is set, all steps required for the "+
                 "computation of the boundary description through f (STEPs 2, 3, and 4) are "+
                 "skipped, and the required parameters for defining f and the scale factor s "+
                 "are read from the file boundary.lis instead, which is always written after "+
                 "the boundary distance from the solute has been adjusted (STEP 4). The reason "+
                 "for this option is that STEP 4 is quite time consuming." ),
        _( "s", type="bool", default=False,
            help="If the -s option is set, the file surface_stat is written, "+
                 "which contains a set of grid points specifying the solvent surface, "+
                 "the error statistics for the distance estimate described above, and "+
                 "information about how many water molecules belong to which group of molecules." ),
        _( "v", type="bool", default=False,
            help="If the -v option is set, the file volume_stat is written (attention: "+
                 "this file may become quite large!), which lists for every grid point "+
                 "within the solvent volume its x-, y-, and z-coordinate, the value f(x,y,z) of "+
                 "the density function at this point, its accurate distance and the "+
                 "approximate distance (which is the efficient estimate used in subsequent "+
                 "MD-simulations) from the solute surface, as well as an approximate value "+
                 "of the curvature of the surface at the point next to (x,y,z)(which also is an efficient "+
                 "estimate that can be used in MD-simulations)." ),
        _( "bulk", type="bool", default=False,
            help="The -bulk option suppresses output of buried water molecules, i.e., "+
                 "only bulk water is written to the output pdb-file." ),
        _( "w", type="bool", default=False,
            help="The -w option suppresses output of the solute, i.e., only water molecules "+
                 "(usually with a hole in the middle, where the solute is located) and ions are "+
                 "written to the output pdb-file." ),
        _( "ion", type="bool", default=False,
            help="If the -ion option is set, STEP 9 is carried out to place sodium and chloride "+
                 "ions into the solvent according to a isotonic Debye-Hueckel density. If the option is "+
                 "not given, STEP 9 is skipped, and no ions are output. psf-file is needed in the same directory as the input." ),
        _( "charge|q", type="float", default=0,
            help="Use this option in addition to -ion to control the total charge of the output system. "+
                 "-q 0 will produce a neutral system." ),
        _( "psf", type="bool", default=False,
            help="If the -psf option is set, SOLVATE writes an X-PLOR-script 'mkpsf.inp' "+
                 "which can be used to generate a structure file for the solute/solvent-system, "+
                 "as required for subsequent MD-simulations (the command xplor < mkpsf.inp will do the job)." ),
    ]
    out = [
        _( "solvate_file", file="{pdb_input.stem}_solvate.pdb" )
    ]
    tmpl_dir = TMPL_DIR

    def _init( self, *args, **kwargs ):
        if self.ion:
            pass
            # TODO: create xplor-file and put it in the same directory as the input
            # equilibration of the output structure can be done with programs,
            # e.g., CHARMm, X-PLOR, or EGO
        extra=[]
        if self.ug: extra.append("-ug")
        if self.ub: extra.append("-ub")
        if self.s: extra.append("-s")
        if self.v: extra.append("-v")
        if self.bulk: extra.append("-bulk")
        if self.w: extra.append("-w")
        if self.ion: extra.append("-ion")
        if self.ion and self.charge:
            extra.append("-q")
            extra.append(self.charge)
        if self.psf: extra.append("-psf")
        self.cmd = [ 
            solvate_CMD,
            "-t", self.thick,
            "-r", self.radius,
            "-n", self.ngauss, 
            self.pdb_input[:-4],
            self.solvate_file[:-4]

        ]
        self.cmd=self.cmd + extra