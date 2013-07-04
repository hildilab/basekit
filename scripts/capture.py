#! /usr/bin/env python

"""A wrapper around the 'capture' web service."""

from basekit.utils.tool import parse_args
from basekit.capture import Capture



def main():
    args, kwargs = parse_args( Capture )
    print Capture( *args, **kwargs )


if __name__ == "__main__":
    main()

