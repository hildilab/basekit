#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.msms import Msms



def main():
    args, kwargs = parse_args( Msms )
    msms = Msms( *args, **kwargs )
    print msms


if __name__ == "__main__":
    main()