#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.dssp import Dssp



def main():
    args, kwargs = parse_args( Dssp )
    print Dssp( *args, **kwargs )


if __name__ == "__main__":
    main()

