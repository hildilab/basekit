#! /usr/bin/env python

"""A collection of MPstruc related tools"""

from basekit.utils.tool import parse_subargs
from basekit.mpstruc import (
	MpstrucDownload, MpstrucInfo
)




def main():
    tools = {
        "db": MpstrucDownload,
        "info": MpstrucInfo
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()