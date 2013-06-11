import os
from string import Template

from utils.tool import CmdTool, ScriptMixin
from utils.job import run_command


DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "jmol" )

JAVA_CMD = "java" 

JMOL_PATH = os.environ.get("JMOL_PATH", "")
JMOL_JAR = os.path.join( JMOL_PATH, "JmolData.jar" )



class Jmol( CmdTool, ScriptMixin ):
    args = [
        { "name": "script_file", "type": "file", "ext": "jspt" }
    ]
    script_tmpl = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, *args, **kwargs ):
        # TODO when args is empty check for self.script_tmpl
        #     to get rid of script_file
        if script_file=="__tmpl__":
            self.script_file = self.outpath( self.script_tmpl )
        self.cmd = [
            JAVA_CMD, "-Xmx4096M", "-jar", JMOL_JAR, "-oxdl", "-s", self.script_file
        ]


class JmolImage( Jmol ):
    args = [
        { "name": "jmol_file", "type": "file", "ext": "jmol" },
        { "name": "scale", "type": "slider", "range": [0, 4], 
            "default": 1, "fixed": True },
        { "name": "width", "type": "slider", "range": [0, 2048], "default": 0 },
        { "name": "height", "type": "slider", "range": [0, 2048], "default": 0 },
        { "name": "cartoon_fancy", "type": "checkbox", "default": True }
    ]
    out = [
        { "name": "image_file", "file": "image.jpg" }
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
        { "name": "dat_file", "type": "file", "ext": "dat", "help": "mrc, obj" },
        { "name": "sigma", "type": "slider", "range": [-1, 5], "default": 0, 
            "fixed": True, "help": "level at which the surface will be created" },
        { "name": "cutoff", "type": "slider", "range": [0, 255], "default": 0, 
            "fixed": True, "help": "value at which the data is ignored" },
        { "name": "resolution", "type": "slider", "range": [-1, 10], 
            "default": 3, "fixed": True }
    ]
    out = [
        { "name": "jvxl_file", "file": "dat.jvxl" }
    ]
    script_tmpl = "jvxl.jspt"
    def _init( self, *args, **kwargs ):
        super(JmolJvxl, self)._init( "__tmpl__" )
    def _pre_exec( self ):
        self._make_script_file(
            dat_file=self.dat_file,
            sigma="" if self.sigma in (None, -1) else "SIGMA %0.2f"%self.sigma,
            cutoff="" if self.cutoff==None else "CUTOFF %0.1f"%self.cutoff,
            resolution="" if self.resolution in (None, -1) else "RESOLUTION %0.1f"%self.resolution,
            jvxl_file=self.jvxl_file
        )


class JmolMovie( Jmol ):
    pass

