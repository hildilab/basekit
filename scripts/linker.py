#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.linker import Linker



def main():
    args, kwargs = parse_args( Linker )
    print Linker( *args, **kwargs )


if __name__ == "__main__":
    main()

