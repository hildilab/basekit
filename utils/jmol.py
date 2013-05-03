

import os

from job import run_command




JMOL_PATH = os.environ["JMOL_PATH"]
JMOL_JAR = os.path.join( JMOL_PATH, "Jmol.jar" )
JMOL_DATA_JAR = os.path.join( JMOL_PATH, "Jmol.jar" )


def run_jmol_script( script, log="jmol.log" ):
	script_file = "jmol.jspt"
	with open( script_file, "w" ) as fp:
		fp.write( script )
	script += "\nexitJmol;"
	cmd = "java -Xmx4096M -jar %s -ionx -s %s" % (
		JMOL_JAR, script_file
	)
	# cmd = "echo -e '%s' | jmol -Ion " % (
	# 	script
	# )
	run_command( cmd, close_fds=True, stdout=False, log=log )


# run_jmol_script( Template( tpl_str ).substitute( **values_dict ) )