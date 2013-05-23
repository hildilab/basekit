#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.jmol import Jmol, JmolImage




def main():
    tools = {
        "script": Jmol, 
        "image": JmolImage
    }
    Tool, args, kwargs = parse_subargs( tools )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()