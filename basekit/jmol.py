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
    jspt_tmpl_file = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, **kwargs ):
        if script_file=="__tmpl__":
            script_file = self.outpath( self.tmpl_file )
        self.script_file = self.abspath( script_file )
        self.cmd = [
            JAVA_CMD, "-Xmx4096M", "-jar", JMOL_JAR, "-oxdl", "-s", self.script_file
        ]


class JmolImage( Jmol ):
    args = [
        { "name": "jmol_file", "type": "file", "ext": "jmol" },
        { "name": "scale", "type": "slider", "range": [0, 4], "default_value": 1, "fixed": True },
        { "name": "width", "type": "slider", "range": [0, 2048], "default_value": 0 },
        { "name": "height", "type": "slider", "range": [0, 2048], "default_value": 0 },
        { "name": "cartoon_fancy", "type": "checkbox", "default_value": True }
    ]
    tmpl_file = "image.jspt"
    def _init( self, jmol_file, scale="", width="", height="", cartoon_fancy=True, **kwargs ):
        self.jmol_file = self.abspath( jmol_file )
        self.width = str(width) if width else "0"
        self.height = str(height) if height else "0"
        self.scale = str(scale) if scale else "0"
        self.cartoon_fancy = "true" if cartoon_fancy else "false"
        self.image_file = self.outpath( "image.jpg" )
        super(JmolImage, self)._init( "__tmpl__" )
        self.output_files = [ "image.jpg" ]
    def _pre_exec( self ):
        self._make_script_file(
            jmol_file=self.jmol_file,
            scale=self.scale,
            width=self.width,
            height=self.height,
            cartoon_fancy=self.cartoon_fancy,
            image_file=self.image_file
        )


class JmolJvxl( Jmol ):
    """Create a JVXL file (the Jmol surface format) from various file formats"""
    args = [
        { "name": "dat_file", "type": "file", "ext": "dat", "help": "mrc, obj" },
        { "name": "sigma", "type": "slider", "range": [0, 5], "default_value": 0, 
            "fixed": True, "help": "level at which the surface will be created" },
        { "name": "cutoff", "type": "slider", "range": [0, 255], "default_value": 0, 
            "fixed": True, "help": "value at which the data is ignored" },
        { "name": "resolution", "type": "slider", "range": [0, 10], "default_value": 3 }
    ]
    tmpl_file = "jvxl.jspt"
    def _init( self, dat_file, sigma=None, cutoff=None, resolution=None, **kwargs ):
        self.dat_file = self.abspath( dat_file )
        self.jvxl_file = self.outpath( "dat.jvxl" )
        self.sigma = sigma
        self.cutoff = cutoff
        self.resolution = resolution
        super(JmolJvxl, self)._init( "__tmpl__" )
        self.output_files = [ self.jvxl_file ]
    def _pre_exec( self ):
        self._make_script_file(
            dat_file=self.dat_file,
            sigma="" if self.sigma==None else "SIGMA %0.2f"%self.sigma,
            cutoff="" if self.cutoff==None else "CUTOFF %0.3f"%self.cutoff,
            resolution="" if self.resolution==None else "RESOLUTION %i"%self.resolution,
            jvxl_file=self.jvxl_file
        )


class JmolMovie( Jmol ):
    pass

