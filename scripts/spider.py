#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.spider import (
    Spider, SpiderConvert, SpiderDeleteFilledDensities, 
    SpiderBox, SpiderCrosscorrelation, LoopCrosscorrel
)




def main():
    tools = {
        "script": Spider, 
        "convert": SpiderConvert, 
        "delFilledDens": SpiderDeleteFilledDensities, 
        "box": SpiderBox, 
        "crosscorrel": SpiderCrosscorrelation,
        "loopcorrel": LoopCrosscorrel
    }
    Tool, args, kwargs = parse_subargs( tools )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()