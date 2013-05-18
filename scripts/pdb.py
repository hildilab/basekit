#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.pdb import PdbDownload, PdbSplit




def main():
    tools = {
        "get": PdbDownload, 
        "split": PdbSplit
    }
    Tool, args, kwargs = parse_subargs( tools )
    tool = Tool( *args, **kwargs )
    print tool



if __name__ == "__main__":
    main()