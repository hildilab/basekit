#! /usr/bin/env python

"""Mapman related tools."""

from basekit.utils.tool import parse_args
from basekit.mapman import Mapman



def main():
    args, kwargs = parse_args( Mapman )
    print Mapman( *args, **kwargs )


if __name__ == "__main__":
    main()