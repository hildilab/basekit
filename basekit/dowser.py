"""
    http://danger.med.unc.edu/hermans/dowser/Dowman.htm
    http://danger.med.unc.edu/hermans/dowser/Cameron_Mura_DOWSER.html
    http://www.ks.uiuc.edu/Research/vmd/plugins/dowser/
    http://muralab.org/~cmura/DOWSER/
    http://www.ks.uiuc.edu/Research/vmd/plugins/dowser/dowser-rna.html
"""

import os
import shutil

from utils import copy_dict
from utils.tool import _, _dir_init, CmdTool, PyTool, ProviMixin
from utils.job import run_command

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "dowser" )

DOWSER_CMD = "dowser"
DOWSERX_CMD = "dowserx"
DOWSER_REPEAT_CMD = "dowser-repeat"

DOWSER_ATOMDICT = os.path.join( TMPL_DIR, "atomdict.db" )
DOWSER_ATOMPARMS = os.path.join( TMPL_DIR, "atomparms.db" )

DOWSER_ARGS = [
    _( "pdb_file", type="file", ext="pdb" ),
    _( "hetero", type="checkbox", default=False ),
    _( "noxtalwater", type="checkbox", default=False ),
    _( "onlyxtalwater", type="checkbox", default=False ),
    _( "probe", type="float", default=0.2 ),
    _( "separation", type="float", default=1.0 ),
    _( "atomtypes", type="text", default=DOWSER_ATOMDICT ),
    _( "atomparms", type="text", default=DOWSER_ATOMPARMS )
]
DOWSER_OUT = [
    _( "input_file", file="myinput.pdb" ),
    _( "wat_file", file="dowserwat.pdb" ),
    _( "watall_file", file="dowserwat_all.pdb" ),
    _( "intsurf_file", file="intsurf.pdb" )
]

class Dowser( CmdTool, ProviMixin ):
    """ A wrapper around the 'dowser' programm. """
    args = DOWSER_ARGS + [
        _( "alt", type="select", options=["x", "repeat"], default=None )
    ]
    out = DOWSER_OUT
    tmpl_dir = TMPL_DIR
    provi_tmpl = "dowser.provi"
    def _init( self, *args, **kwargs ):
        exe = DOWSER_CMD
        if self.alt=="x":
            exe = DOWSERX_CMD
        elif self.alt=="repeat":
            exe = DOWSER_REPEAT_CMD

        self.cmd = [ 
            exe, self.relpath( self.input_file ), 
            "-probe", self.probe,
            "-separation", self.separation
        ]
        if self.hetero: self.cmd.append( "-hetero" )
        if self.noxtalwater: self.cmd.append( "-noxtalwater" )
        if self.onlyxtalwater: self.cmd.append( "-onlyxtalwater" )
        if self.atomtypes: self.cmd += [ "-atomtypes", self.atomtypes ]
        if self.atomparms: self.cmd += [ "-atomparms", self.atomparms ]
    def _pre_exec( self ):
        shutil.copy( self.pdb_file, self.input_file )
    def _post_exec( self ):
        self._make_provi_file(
            input_file=self.relpath( self.input_file ),
            wat_file=self.relpath( self.wat_file ),
            watall_file=self.relpath( self.watall_file ),
            intsurf_file=self.relpath( self.intsurf_file )
        )


class DowserRepeat( PyTool, ProviMixin ):
    """ A wrapper around the 'dowser-repeat' programm, which is
        called until no more new waters are found.
    """
    args = DOWSER_ARGS + [
        _( "alt", type="select", options=["x"], default=None )
    ]
    out = DOWSER_OUT + [
        _( "repeat_dir", dir="repeats" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "dowser.provi"
    def _init( self, *args, **kwargs ):
        self.kwargs = kwargs
    def func( self ):
        """ Initial procedure by Dominic Theune. """
        rep_count = 0
        dowserwat = []
        dowserwat_all = []
        while True:
            alt = "repeat" if rep_count else self.alt
            rep_out = os.path.join( self.repeat_dir, str(rep_count) )
            dowser = Dowser( self.pdb_file, **copy_dict( 
                self.kwargs, run=False, output_dir=rep_out, alt=alt
            ))
            # write dowser waters found until now
            with open( dowser.wat_file, "w" ) as fp:
                fp.writelines( dowserwat )
            with open( dowser.watall_file, "w" ) as fp:
                fp.writelines( dowserwat_all )
            # exec dowser
            dowser()
            # read newly found dowser waters
            with open( dowser.wat_file, "r" ) as fp:
                new_wat = fp.readlines()
            with open( dowser.watall_file, "r" ) as fp:
                new_watall = fp.readlines()
            # checks if there are new waters
            if len( new_wat ) > 1 and len( new_watall ) > 1:
                dowserwat += new_wat
                dowserwat_all += new_watall
                rep_count += 1
            else:
                break
        # write every dowser water found
        with open( self.wat_file, "w" ) as fp:
            fp.writelines( dowserwat )
        with open( self.watall_file, "w" ) as fp:
            fp.writelines( dowserwat_all )
        # copy intsurf
        shutil.copy( dowser.intsurf_file, self.intsurf_file )
    def _pre_exec( self ):
        shutil.copy( self.pdb_file, self.input_file )
    def _post_exec( self ):
        self._make_provi_file(
            input_file=self.relpath( self.input_file ),
            wat_file=self.relpath( self.wat_file ),
            watall_file=self.relpath( self.watall_file ),
            intsurf_file=self.relpath( self.intsurf_file )
        )

