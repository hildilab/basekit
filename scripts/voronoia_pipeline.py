#! /usr/bin/env python

"""Voronoia pipeline."""

from basekit.utils.tool import parse_subargs
from basekit.voronoia_pipeline import VoronoiaPipeline



def main():
    tools = {
        "pipe": VoronoiaPipeline
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()