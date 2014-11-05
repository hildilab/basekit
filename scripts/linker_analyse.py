#! /usr/bin/env python

"""A collection of tools loops linking residues."""

from basekit.utils.tool import parse_subargs
from basekit.linker_analyse import(
    LinkIt_analyse, LinkItDensity, LnkItVali, AnalyseLiniktRun,CutPDB,Cut_analyse_PDB2,LinkItDensity2_analyse,CutPDB3,LinkItDensity3
    #MultiLinkItLinkerTest, 
)



def main():
    tools = {
        "linkit": LinkIt_analyse, 
        "linkit+dens": LinkItDensity,
        #"test": LinkerTest,
        "linkvali": LnkItVali,
        "analyse" : AnalyseLiniktRun,
        #"multi-linkit": MultiLinkIt
        "cut":CutPDB,
        "cut2_analyse":Cut_analyse_PDB2,
        "linkit+dens2":LinkItDensity2_analyse,
        "cut3":CutPDB3,
        "linkit+dens3":LinkItDensity3,        
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()

