#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.motif import CapsMotifFinder, CapsMotifPlotter



def main():
    tools = {
        "finder": CapsMotifFinder,
        "plotting": CapsMotifPlotter
    }
    Tool, args, kwargs = parse_subargs( tools )
    Tool( *args, **kwargs )



if __name__ == "__main__":
    main()

