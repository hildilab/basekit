#! /usr/bin/env python

"""A collection of wrappers around the tmdet program."""

from basekit.utils.tool import parse_subargs
from basekit.tmdet import (
	Tmdet, PdbtmInfo, Pdbtm, PdbtmList, PdbtmDownload
)


def main():
    tools = {
        "pdb": Tmdet,
        "info": PdbtmInfo,
        "id": Pdbtm,
        "list": PdbtmList,
        "db": PdbtmDownload
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()