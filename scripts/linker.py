#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.linker import LinkIt



def main():
    args, kwargs = parse_args( LinkIt )
    print LinkIt( *args, **kwargs )


if __name__ == "__main__":
    main()

