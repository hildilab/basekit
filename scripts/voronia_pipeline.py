#! /usr/bin/env python

"""Voronoia pipeline."""

from basekit.utils.tool import parse_args
from basekit.voronoia_pipeline import Voronoia



def main():
    tools = {
        "finder": CapsMotifFinder,
        "strucinfo": StructureGetter,
        "superpose": Superpose
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()