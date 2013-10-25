#! /usr/bin/env python

"""A collection of pdb related tools."""

from basekit.utils.tool import parse_subargs
from basekit.pdb import (
    PdbDownload, PdbSplit, PdbUnzip, PdbHetDictionary, 
    PdbSuperpose, PdbEdit, PdbInfo, NumpdbTest, RnaList,
    ListCompare, ListJoin, PdbAssembly,SplitPdbSSE, LoopDelete,
    MakeAllRotameres
)




def main():
    tools = {
        "SplitSSE": SplitPdbSSE,
        "get": PdbDownload, 
        "split": PdbSplit,
        "unzip": PdbUnzip,
        "het": PdbHetDictionary,
        "superpose": PdbSuperpose,
        "edit": PdbEdit,
        "info": PdbInfo,
        "test": NumpdbTest,
        "rna": RnaList,
        "compare": ListCompare,
        "join": ListJoin,
        "assembly": PdbAssembly,
        "delete": LoopDelete,
	"rotamere": MakeAllRotameres
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()