#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.msms import Msms



def main():
    args, kwargs = parse_args( Msms )
    Msms( *args, **kwargs )


if __name__ == "__main__":
    main()