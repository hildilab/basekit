#! /usr/bin/env python

"""Mapman related tools."""

from basekit.utils.tool import parse_args, parse_subargs
from basekit.mapman import Mapman, BrixToMap



def main():
    tools = {
        "mapman": Mapman,
        "b2m": BrixToMap
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )

if __name__ == "__main__":
    main()