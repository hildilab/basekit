import os
from string import Template

from utils.tool import _, _dir_init, CmdTool, ScriptMixin
from utils.job import run_command

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "jmol" )

JAVA_CMD = "java" 

JMOL_PATH = os.environ.get("JMOL_PATH", "")
JMOL_JAR = os.path.join( JMOL_PATH, "JmolData.jar" )



class Jmol( CmdTool, ScriptMixin ):
    args = [
        _( "script_file", type="file", ext="jspt" )
    ]
    script_tmpl = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, *args, **kwargs ):
        # TODO when args is empty check for self.script_tmpl
        #     to get rid of script_file
        if script_file=="__tmpl__":
            self.script_file = self.outpath( self.script_tmpl )
        self.cmd = [
            JAVA_CMD, "-Xmx6000M", "-jar", JMOL_JAR, "-oxdl", "-s", self.script_file
        ]


class JmolImage( Jmol ):
    args = [
        _( "jmol_file", type="file", ext="jmol" ),
        _( "scale", type="float", range=[0, 8], step=1, default=2 ),
        _( "width", type="int", range=[0, 4096], default=0 ),
        _( "height", type="int", range=[0, 4096], default=0 ),
        _( "cartoon_fancy", type="bool", default=True )
    ]
    out = [
        _( "image_file", file="image.jpg" )
    ]
    script_tmpl = "image.jspt"
    def _init( self, *args, **kwargs ):
        super(JmolImage, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file(
            jmol_file=self.jmol_file,
            scale=self.scale or "0",
            width=self.width or "0",
            height=self.height or "0",
            cartoon_fancy="true" if self.cartoon_fancy else "false",
            image_file=self.image_file
        )


class JmolJvxl( Jmol ):
    """Create a JVXL file (the Jmol surface format) from various file formats"""
    args = [
        _( "dat_file", type="file", ext="dat", help="mrc, obj, mrc, ccp4" ),
        _( "struc_file", type="file", ext="pdb", default=None ),
        _( "sigma", type="slider", range=[-1, 5], default=0, 
            fixed=True, help="level at which the surface will be created" ),
        _( "within", type="slider", range=[-1, 5], default=None, 
            fixed=True, help="within distance of atoms" ),
        _( "cutoff", type="slider", range=[0, 255], default=None, 
            fixed=True, help="value at which the data is ignored" ),
        _( "resolution", type="slider", range=[-1, 10], 
            default=3, fixed=True ),
        _( "style", type="str", default=None )
    ]
    out = [
        _( "jvxl_file", file="dat.jvxl" )
    ]
    script_tmpl = "jvxl.jspt"
    def _init( self, *args, **kwargs ):
        super(JmolJvxl, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file(
            dat_file=self.dat_file,
            struc_file="" if self.struc_file==None else self.struc_file,
            sigma="" if self.sigma in (None, -1) else "SIGMA %0.2f"%self.sigma,
            within="" if self.within==None else "WITHIN %0.2f {*}"%self.within,
            cutoff="" if self.cutoff==None else "CUTOFF %0.1f"%self.cutoff,
            resolution="" if self.resolution in (None, -1) else "RESOLUTION %0.1f"%self.resolution,
            style="" if self.style==None else self.style,
            jvxl_file=self.jvxl_file
        )


class JmolMovie( Jmol ):
    pass

