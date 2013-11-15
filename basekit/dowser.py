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
from utils.numpdb import NumPdb


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "dowser" )

DOWSER_CMD = "dowser"
DOWSERX_CMD = "dowserx"
DOWSER_REPEAT_CMD = "dowser-repeat"

DOWSER_ATOMDICT = os.path.join( TMPL_DIR, "atomdict.db" )
DOWSER_ATOMPARMS = os.path.join( TMPL_DIR, "atomparms.db" )


class DowserMixin( object ):
    def _pre_exec( self ):
        shutil.copy( self.pdb_file, self.input_file )
    def _post_exec( self ):
        self._make_provi_file(
            input_file=self.relpath( self.input_file ),
            wat_file=self.relpath( self.wat_file ),
            watall_file=self.relpath( self.watall_file ),
            intsurf_file=self.relpath( self.intsurf_file )
        )
        # make sure all file exist
        for f in [ self.wat_file, self.watall_file, self.intsurf_file ]:
            with open( f, "a" ):
                pass
        # write a pdb with all dowser waters but no others
        with open( self.dowser_file, "w" ) as fp:
            with open( self.pdb_file, "r" ) as fp_pdb:
                for line in fp_pdb:
                    if ( line[0:6] in ["ATOM  ", "HETATM"] and
                            line[17:20]!="HOH" ):
                        fp.write( line )
            with open( self.wat_file, "r" ) as fp_wat:
                for line in fp_wat:
                    if line[0:6] in ["HETATM"]:
                        fp.write( line )
        # rename HOH.OW atoms to HOH.O and
        # repair atom index
        npdb = NumPdb( self.dowser_file )
        for i, a in enumerate( npdb._atoms, start=1 ):
            if a["atomname"]==" OW " and a["resname"]=="HOH":
                a["atomname"] = " O  "
            a["atomno"] = i
        npdb.write( self.dowser_file )


DOWSER_ARGS = [
    _( "pdb_file", type="file", ext="pdb" ),
    _( "hetero", type="bool", default=False ),
    _( "noxtalwater", type="bool", default=False ),
    _( "onlyxtalwater", type="bool", default=False ),
    _( "probe", type="float", default=0.4 ),
    _( "separation", type="float", default=1.0 ),
    _( "atomtypes", type="str", default=DOWSER_ATOMDICT ),
    _( "atomparms", type="str", default=DOWSER_ATOMPARMS )
]
DOWSER_OUT = [
    _( "input_file", file="myinput.pdb" ),
    _( "wat_file", file="dowserwat.pdb" ),
    _( "watall_file", file="dowserwat_all.pdb" ),
    _( "intsurf_file", file="intsurf.pdb" ),
    _( "dowser_file", file="dowser.pdb" )
]

class Dowser( DowserMixin, CmdTool, ProviMixin ):
    """ A wrapper around the 'dowser' programm. """
    args = DOWSER_ARGS + [
        _( "alt", type="str", options=["x", "repeat"], default=None )
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
    

class DowserRepeat( DowserMixin, PyTool, ProviMixin ):
    """ A wrapper around the 'dowser-repeat' programm, which is
        called until no more new waters are found.
    """
    args = DOWSER_ARGS + [
        _( "alt", type="str", options=["x"], default=None ),
        _( "max_repeats", type="int", default=None )
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
        max_resno = 0
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
            # increment water resno
            new_wat2 = []
            for l in new_wat:
                if l[0:6]=="HETATM":
                    resno = int( l[22:26] ) + max_resno
                    l = l[0:22] + "{:>4}".format( resno ) + l[26:]
                new_wat2.append( l )
            new_watall2 = []
            for l in new_watall:
                if l[0:6]=="HETATM":
                    resno = int( l[22:26] ) + max_resno
                    l = l[0:22] + "{:>4}".format( resno ) + l[26:]
                new_watall2.append( l )
            if new_watall2:
                max_resno = max( map( lambda x: int(x[22:26]), new_watall2 ) )
            # check if there are new waters
            if ( len( new_wat2 )>1 and len( new_watall2 )>1 and
                    ( not self.max_repeats or rep_count<self.max_repeats-1 ) ):
                dowserwat += new_wat2
                dowserwat_all += new_watall2
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


