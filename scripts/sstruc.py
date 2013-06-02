#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.sstruc import Sstruc, SstrucFinder, SstrucTest, SstrucParallel




def main():
    tools = {
    	"test": SstrucTest,
        "pdb": Sstruc,
        "find": SstrucFinder,
        "parallel": SstrucParallel
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()