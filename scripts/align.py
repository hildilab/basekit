#! /usr/bin/env python

"""Alignment tools."""

from basekit.utils.tool import parse_args
from basekit.align import Muscle



def main():
    args, kwargs = parse_args( Muscle )
    print Muscle( *args, **kwargs )


if __name__ == "__main__":
    main()