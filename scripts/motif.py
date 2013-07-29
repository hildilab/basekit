#! /usr/bin/env python


"""A collection of tools to find cap motifs and plot their dihedral angles. """


from basekit.utils.tool import parse_subargs
from basekit.motif import CapsMotifFinder, StructureGetter, Superpose



def main():
    tools = {
        "finder": CapsMotifFinder,
        "strucinfo": StructureGetter,
        "superpose": Superpose
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()

