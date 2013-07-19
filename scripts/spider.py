#! /usr/bin/env python

"""A collection of tools based on spider scripts."""

from basekit.utils.tool import parse_subargs
from basekit.spider import (
    MrcHeaderPrinter,
    Spider, SpiderShift, SpiderConvert, SpiderDeleteFilledDensities, 
    SpiderBox, SpiderReConvert, SpiderCrosscorrelation, LoopCrosscorrel, 
    SpiderPdbBox
)



def main():
    tools = {
        "mrc": MrcHeaderPrinter,
        "script": Spider,
        "shift": SpiderShift,
        "convert": SpiderConvert, 
        "delFilledDens": SpiderDeleteFilledDensities, 
        "box": SpiderBox,
        "reconvert": SpiderReConvert,
        "crosscorrel": SpiderCrosscorrelation,
        "loopcorrel": LoopCrosscorrel,
        "pdbbox": SpiderPdbBox
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()