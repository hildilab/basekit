#! /usr/bin/env python

"""Voronoia related tools."""

from basekit.utils.tool import parse_args
from basekit.voronoia import Voronoia



def main():
    args, kwargs = parse_args( Voronoia )
    print Voronoia( *args, **kwargs )


if __name__ == "__main__":
    main()

