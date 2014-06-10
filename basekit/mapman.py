from __future__ import with_statement
from __future__ import division



import os
import logging

from basekit import utils
from utils import path
from utils.tool import _, _dir_init, CmdTool

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "mapman" )
Mapman_CMD = os.path.join( TMPL_DIR, "lx_mapman" )

logging.basicConfig()
LOG = logging.getLogger('mapman')
LOG.setLevel( logging.ERROR )



class Mapman( CmdTool ):
    """A wrapper around the 'Mapman' programm."""
    args = [
        _( "em_input", type="file", help="Input EM-File of type [PROTEIN FFT-Y TENEYCK2 CCP4 X-PLOR OLDEZD MASK NEWEZD BINXPLOR BRICK DSN6 3DMATRIX TNT PHASES FSMASK BRIX XPLOR CNS EZD EM08 OMAP MPI AMBER]" ),
        _( "newformat|nf", type="str", default="map", 
            options=['ccp4', 'oldezd', 'mask', 'newezd', 'envelope',
                     'x-plor', 'dsn6', 'brix', 'xplor', 'cns', 'ezd',
                     'amore', 'omap', 'turbo', 'mpi'],
            help="output format of type [CCP4 OLDEZD MASK NEWEZD ENVELOPE X-PLOR DSN6 BRIX XPLOR CNS EZD AMORE OMAP TURBO MPI]" ),

    ]
    out = [
        _( "newbase", file="{em_input.stem}.{newformat}" ),
        _( "params_file", file="{em_input.stem}_mapman.params" ),

    ]
    tmpl_dir = TMPL_DIR
    def _init( self, *args, **kwargs ):
        self.newformat = self.newformat.upper()
        if self.newformat=='MAP':
            self.newformat='CCP4'
        self.format=utils.path.ext(self.em_input)[1:].upper()
        if self.format=='MAP':
            self.format='CCP4'
        self.cmd = [ Mapman_CMD ]
        self.cmd_input = self.params_file
        self.no_cmd = False
    def _pre_exec( self ):
        add=""
        if self.newformat=='MASK':
            add="0.0010"
        with open(self.params_file, "w") as fp:
            fp.write( "zp 90000000 2\nre slot "+\
                self.em_input+\
                " "+\
                self.format+\
                " \nwr slot "+\
                self.relpath( self.newbase )+\
                " "+\
                self.newformat+\
                add+\
                " \nquit\n"
            )
    def _post_exec( self ):
        if self.check( full=True ):
            pass
                

