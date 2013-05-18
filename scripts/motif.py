#! /usr/bin/env python

from basekit.utils.tool import parse_args
from basekit.motif import CapsMotifFinder



def main():
    args, kwargs = parse_args( CapsMotifFinder )
    caps_motif_finder = CapsMotifFinder( *args, **kwargs )
    print caps_motif_finder


if __name__ == "__main__":
    main()

