#! /usr/bin/env python

"""A collection of OPM related tools"""

from basekit.utils.tool import parse_subargs
from basekit.opm import (
	Opm, OpmList, OpmInfo, Ppm, Ppm2
)




def main():
    tools = {
        "id": Opm,
        "list": OpmList,
        "info": OpmInfo,
        "pdb": Ppm,
        "id2": Ppm2
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()