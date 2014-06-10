#! /usr/bin/env python

"""Alignment tools."""

from basekit.utils.tool import parse_args, parse_subargs
from basekit.align import Muscle, TheseusMakeFasta



def main():
    #args, kwargs = parse_args( Muscle )
    #print Muscle( *args, **kwargs )
    tools = {
        "muscle": Muscle,
        "tmf": TheseusMakeFasta
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )

if __name__ == "__main__":
    main()