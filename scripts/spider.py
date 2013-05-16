#! /usr/bin/env python

from basekit.utils.tool import parse_subargs
from basekit.spider import Spider, SpiderConvert, SpiderDeleteFilledDensities, SpiderBox, SpiderCrosscorrelation




def main():
    tools = {
        "script": Spider, 
        "convert": SpiderConvert, 
        "delFilledDens": SpiderDeleteFilledDensities, 
        "box": SpiderBox, 
        "crosscorrel": SpiderCrosscorrelation
    }
    Tool, args, kwargs = parse_subargs( tools )
    Tool( *args, **kwargs )



if __name__ == "__main__":
    main()