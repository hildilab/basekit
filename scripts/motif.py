#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.motif import CapsMotifFinder



def main():
    args, kwargs = parse_args( CapsMotifFinder )
    CapsMotifFinder( *args, **kwargs )


if __name__ == "__main__":
    main()

