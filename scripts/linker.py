#! /usr/bin/env python

"""A collection of tools loops linking residues."""

from basekit.utils.tool import parse_subargs
from basekit.linker import(
    LinkIt, LinkItDensity, LnkItVali, AnalyseLiniktRun,CutPDB
    #MultiLinkItLinkerTest, 
)



def main():
    tools = {
        "linkit": LinkIt, 
        "linkit+dens": LinkItDensity,
        #"test": LinkerTest,
        "linkvali": LnkItVali,
        "analyse" : AnalyseLiniktRun,
        #"multi-linkit": MultiLinkIt
        "cut":CutPDB
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()

