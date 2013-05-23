#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.apbs import Apbs



def main():
    args, kwargs = parse_args( Apbs )
    print Apbs( *args, **kwargs )


if __name__ == "__main__":
    main()

