#! /usr/bin/env python

"""A collection of pdb related tools."""

from basekit.utils.tool import parse_subargs
from basekit.pdb import (
	PdbDownload, PdbSplit, PdbUnzip, PdbHetDictionary, 
	PdbSuperpose, PdbBox, NumpdbTest
)




def main():
    tools = {
        "get": PdbDownload, 
        "split": PdbSplit,
        "unzip": PdbUnzip,
        "het": PdbHetDictionary,
        "superpose": PdbSuperpose,
        "box": PdbBox,
        "test": NumpdbTest
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()