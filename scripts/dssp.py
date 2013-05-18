#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.dssp import Dssp



def main():
    args, kwargs = parse_args( Dssp )
    dssp = Dssp( *args, **kwargs )
    print dssp


if __name__ == "__main__":
    main()

