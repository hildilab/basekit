#! /usr/bin/env python

"""A collection of tools loops linking residues."""

from basekit.utils.tool import parse_subargs
from basekit.linker_analyse import(
    LinkIt, LinkItDensity, LnkItVali, CutPDB2,AnalyseLiniktRun,CutPDB,Cut_analyse_PDB2,LinkItDensity2_analyse,
    CutPDB3,LinkItDensity3,CalcSheets,MultiLinkItDensity,MultiCutPDB,LinkItList,PdbSequence
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
        "linkit+dens2":LinkItDensity2_analyse,
        "cut3":CutPDB3,
        "linkit+dens3":LinkItDensity3,
        "calcsheets":CalcSheets,
        "multilinkdens":MultiLinkItDensity,
        "multicut":MultiCutPDB,
        "linkitlist":LinkItList,
        "seqlist":PdbSequence
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )


if __name__ == "__main__":
    main()

