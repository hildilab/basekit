#! /usr/bin/env python

"""A collection of wrappers around the dowser program."""

from basekit.utils.tool import parse_subargs
from basekit.dowser import Dowser, DowserRepeat




def main():
    tools = {
        "pdb": Dowser, 
        "repeat": DowserRepeat
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()