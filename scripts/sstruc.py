#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.sstruc import Sstruc, SstrucFinder




def main():
    tools = {
        "info": Sstruc, 
        "find": SstrucFinder
    }
    Tool, args, kwargs = parse_subargs( tools )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()