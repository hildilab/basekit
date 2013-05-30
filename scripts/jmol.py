#! /usr/bin/env python

"""A collection of tools based on Jmol scripts"""

from basekit.utils.tool import parse_subargs
from basekit.jmol import Jmol, JmolImage, JmolJvxl




def main():
    tools = {
        "script": Jmol, 
        "image": JmolImage,
        "jvxl": JmolJvxl
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()