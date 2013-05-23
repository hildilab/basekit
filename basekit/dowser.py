#! /usr/bin/env python

import argparse
import shutil

from utils import dir_walker, working_directory
from utils.job import run_command


DOWSER_CMD = "dowser"
DOWSERX_CMD = "dowserx"
DOWSER_REPEAT_CMD = "dowser-repeat"


from utils.tool import CmdTool


class Dowser( CmdTool ):
    pass


class DowserRepeat( CmdTool ):
    pass



def dowser( input_file, output_dir, log=None,
            hetero=False, noxtalwater=False, onlyxtalwater=False,
            atomtypes=None, atomparms=None, probe=0.2, separation=1.0,
            alt=None ):
    exe = DOWSER_CMD
    if alt=="x":
        exe = DOWSERX_CMD
    elif alt=="repeat":
        exe = DOWSER_REPEAT_CMD

    shutil.copy( input_file, Path( output_dir, "myinput.pdb" ) )
    cmd = "%s %s -probe %f -separation %f" % (
        exe, "myinput.pdb", probe, separation
    )
    if hetero: cmd += " -hetero"
    if noxtalwater: cmd += " -noxtalwater"
    if onlyxtalwater: cmd += " -onlyxtalwater"
    if atomtypes: cmd += " -atomtypes %s" % atomtypes
    if atomparms: cmd += " -atomparms %s" % atomparms
    with working_directory( output_dir ):
        run_command( cmd, close_fds=True, stdout=log )


def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument(
        '-i', '--inputFile', metavar='inputFile',
        help='input file path')
    parser.add_argument(
        '-dw', '--dirWalker', metavar='dirWalker',
        help='dir walker regex pattern e.g. "(.*)_OPM_mod\.pdb"')
    parser.add_argument(
        '-o', '--outputDir', metavar='outputDir', default='.',
        help='the output directory')
    parser.add_argument(
        '-het', '--hetero', metavar='hetero', type=int, default=1,
        help='use hetero atoms')
    parser.add_argument(
        '-sep', '--separation', metavar='separation', type=float, default=1.0,
        help='water separation in angstrom')
    parser.add_argument(
        '-alt', '--alternateUsage', metavar='alternateUsage', type=str, default="",
        help='alternate usage: x for dowserx and repeat for dowser-repeat')

    
    # parse the command line
    args = parser.parse_args()
    
    if args.inputFile:
        input_file = Path( args.inputFile ).expand().absolute()

    # kwargs = {
    #     "hetero": False, "noxtalwater": False, "onlyxtalwater": False,
    #     "atomtypes": None, "atomparms": None, "probe":0.2, "separation": 1.0,
    #     "log": "stdlog.txt"
    # }

    kwargs = {
        "hetero": args.hetero, "noxtalwater": False, "onlyxtalwater": False,
        "atomtypes": None, "atomparms": None, "probe":0.2, "separation": args.separation,
        "log": "stdlog.txt", "alt": args.alternateUsage
    }

    
    if args.inputFile and args.dirWalker:
        print args.dirWalker
        for m, dw_input_file in dir_walker( input_file, args.dirWalker ):
            dw_output_dir = Path( dw_input_file.parent, args.outputDir )
            dw_output_dir.mkdir(True)
            dw_output_dir = dw_output_dir.expand()
            dowser( dw_input_file, dw_output_dir, **kwargs )
    elif args.inputFile:
        output_dir = Path(args.outputDir)
        output_dir.mkdir(True)
        output_dir = output_dir.expand()
        dowser( input_file, output_dir, **kwargs )



if __name__ == "__main__":
    main()



