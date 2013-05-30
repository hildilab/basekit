#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.linker import LinkIt, LinkItDensity, LinkerTest



def main():
    tools = {
        "linkit": LinkIt, 
        "linkit+dens": LinkItDensity,
        "test": LinkerTest
    }
    Tool, args, kwargs = parse_subargs( tools )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()

