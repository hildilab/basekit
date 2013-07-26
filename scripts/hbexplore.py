#! /usr/bin/env python

"""A wrapper around the 'hbexplore' programm."""

from basekit.utils.tool import parse_args
from basekit.hbexplore import HBexplore



def main():
    args, kwargs = parse_args( HBexplore )
    print HBexplore( *args, **kwargs )


if __name__ == "__main__":
    main()

