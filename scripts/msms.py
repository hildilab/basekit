#! /usr/bin/env python

"""A wrapper around the 'msms' programm"""

from basekit.utils.tool import parse_args
from basekit.msms import Msms



def main():
    args, kwargs = parse_args( Msms )
    print Msms( *args, **kwargs )


if __name__ == "__main__":
    main()