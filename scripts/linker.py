#! /usr/bin/env python

"""A collection of tools loops linking residues."""

from basekit.utils.tool import parse_subargs
from basekit.linker import LinkIt, LinkItDensity, LinkerTest, LinkItDensityjason,LnkItVali



def main():
    tools = {
        "linkit": LinkIt, 
        "linkit+dens": LinkItDensity,
        "test": LinkerTest,
        "linkjson": LinkItDensityjason,
        "linkvali": LnkItVali
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()

