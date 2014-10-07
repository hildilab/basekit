#! /usr/bin/env python

"""A collection of tools based on spider scripts."""

from basekit.utils.tool import parse_subargs
from basekit.spider_analyse import (
    MrcHeaderPrinter,
    Spider, SpiderShift, SpiderConvert, SpiderDeleteFilledDensities, 
    SpiderBox, SpiderReConvert, SpiderCrosscorrelation, LoopCrosscorrel, 
    SpiderPdbBox,SpiderSidechainCorrelation,  LoopSidechainCorrelation,SpiderDeleteBackbone,BuildBest,OptimizeRotamer,
    LoopRotamerOptimize,SideChainStatistics, SpiderCropMap,SpiderCropMrc,SpiderAnalyse,LoopCrosscorrel2,LoopCrosscorrel3,SpiderAnalyse2,SpiderAnalyse3,Calcoricc,Looporicc,BuildModel
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
        "pdbbox": SpiderPdbBox,
        "sidechaincc": SpiderSidechainCorrelation,
        "loopsidechain": LoopSidechainCorrelation,
        "deletebb":SpiderDeleteBackbone,
        "buildbest":BuildBest,
        "optimize":OptimizeRotamer,
        "loopoptimize":LoopRotamerOptimize,
        "sidestat":SideChainStatistics,
        "crop":SpiderCropMap,
        "cropmrc":SpiderCropMrc,
        "analyse":SpiderAnalyse,
        "analyse2":SpiderAnalyse2,
        "analyse3":SpiderAnalyse3,
        "loopcorrel2": LoopCrosscorrel2,
        "loopcorrel3": LoopCrosscorrel3,
        "ccori":Calcoricc,
        "loopcc":Looporicc,
        "build":BuildModel
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    print Tool( *args, **kwargs )



if __name__ == "__main__":
    main()