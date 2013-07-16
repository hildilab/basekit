#! /usr/bin/env python

"""A collection of project managment tools."""

from basekit.utils.tool import parse_subargs
from basekit.project import (
    ProjectInfo, ProjectRun
)



def main():
    tools = {
        "info": ProjectInfo,
        "run": ProjectRun,
    }
    Tool, args, kwargs = parse_subargs( tools, description=__doc__ )
    Tool( *args, **kwargs )


if __name__ == "__main__":
    main()

