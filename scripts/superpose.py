#! /usr/bin/env python

"""Theseus related tools."""

from basekit.utils.tool import parse_args
from basekit.superpose import Theseus



def main():
    args, kwargs = parse_args( Theseus )
    print Theseus( *args, **kwargs )


if __name__ == "__main__":
    main()