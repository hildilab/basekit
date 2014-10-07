#! /usr/bin/env python

"""A collection of tools loops linking residues."""

from basekit.utils.tool import parse_subargs
from basekit.linker import(
    LinkIt, LinkItDensity, LnkItVali, AnalyseLiniktRun,CutPDB,CutPDB2,LinkItDensity2,CutPDB3,LinkItDensity3
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
        "cut":CutPDB,
        "cut2":CutPDB2,
        "linkit+dens2":LinkItDensity2,
        "cut3":CutPDB3,
        "linkit+dens3":LinkItDensity3,        
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()

