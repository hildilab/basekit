#! /usr/bin/env python

"""Moderna related tools."""

from basekit.utils.tool import parse_args
from basekit.moderna_tools import Moderna



def main():
    args, kwargs = parse_args( Moderna )
    print Moderna( *args, **kwargs )


if __name__ == "__main__":
    main()
