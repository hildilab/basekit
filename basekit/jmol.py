import os
from string import Template

from utils.tool import CmdTool, ScriptMixin
from utils.job import run_command


DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "jmol" )

JAVA_CMD = "java" 

JMOL_PATH = os.environ.get("JMOL_PATH", "")
JMOL_JAR = os.path.join( JMOL_PATH, "Jmol.jar" )



class Jmol( CmdTool, ScriptMixin ):
    args = [
        { "name": "script_file", "type": "file", "ext": "jspt" }
    ]
    jspt_tmpl_file = None
    tmpl_dir = TMPL_DIR
    def _init( self, script_file, **kwargs ):
        self.script_file = os.path.abspath( script_file )
        self.cmd = [
            JAVA_CMD, "-Xmx4096M", "-jar", JMOL_JAR, "-ionx", "-s", self.script_file
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
    	width = str(width) if width else "0"
    	height = str(height) if height else "0"
    	scale = str(scale) if scale else "0"
    	cartoon_fancy = "true" if cartoon_fancy else "false"
        script_file = self._make_script_file(
        	jmol_file=os.path.abspath( jmol_file ),
        	scale=scale,
        	width=width,
        	height=height,
        	cartoon_fancy=cartoon_fancy
    	)
        super(JmolImage, self)._init( script_file )
        self.output_files = [ "image.jpg" ]


class JmolMovie( Jmol ):
	pass

