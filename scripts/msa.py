#! /usr/bin/env python

"""A collection of multiple sequence alignment related tools."""

from basekit.utils.tool import parse_subargs
from basekit.msa import (
	Muscle
)




def main():
    tools = {
        "muscle": Muscle
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()