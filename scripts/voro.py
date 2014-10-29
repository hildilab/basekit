#! /usr/bin/env python

"""A collection of Voronoia tools"""

from basekit.utils.tool import parse_subargs
from basekit.voro import VoronoiaPipeline, VoronoiaStats




def main():
    tools = {
        "pdb": VoronoiaPipeline,
        "stats": VoronoiaStats
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()