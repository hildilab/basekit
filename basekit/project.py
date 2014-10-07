from __future__ import with_statement

import os
import re
import copy
import importlib
import collections
import json
import tempfile

from utils import working_directory, listify
from utils.tool import _, _dir_init, PyTool


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "project" )


def decode_project_json( pairs ):
    source = None
    key = None
    concat = None
    update = None
    for k, v in pairs:
        if k == '__source__':
            source = v
        elif k == '__key__':
            key = v
        elif k == '__concat__':
            concat = v
        elif k == '__update__':
            update = v
    if source:
        data = read_project_file( source )
        if key:
            data = data[ key ]
        if concat:
            data += concat
        elif update:
            data.update( update )
        return data
    else:
        return collections.OrderedDict( pairs )


def read_project_file( project_file ):
    temp = tempfile.SpooledTemporaryFile()
    with open( project_file, "r" ) as fp:
        for line in fp:
            if not line.strip().startswith("//"):
                temp.write( line )
    temp.seek(0)
    try:
        with working_directory( os.path.dirname( project_file ) ):
            return json.load(
                temp, object_pairs_hook=decode_project_json
            )
    except ValueError as e:
        print "Error reading project json file (%s). %s" % (
            os.path.basename( project_file ), e
        )


def get_tool( t ):
    module, tool = t.get("__name__").rsplit(".", 1)
    return importlib.import_module( module ).__dict__[ tool ]


def get_wd( project_json, p, **kwargs ):
    wd = os.path.join(
        project_json["project"].get("__dir__", ""),
        p["__dir__"]
    )
    return wd.format( **kwargs )


def project_info( project_json ):
    print "{:#^80}".format( " PARTS " )
    for pid, p in project_json["parts"].iteritems():
        print "ID: {id:<25}DIR: {__dir__}".format( id=pid, **p )
        if os.path.isdir( get_wd( project_json, p ) ):
            print " `-- directory available."
        else:
            print " `-- Error: directory not found."

    print "{:#^80}".format( " TOOLS " )
    for tid, t in project_json["tools"].iteritems():
        print "ID: {id:<45}NAME: {__name__}".format( id=tid, **t )
        try:
            get_tool( t )
            print " `-- tool available."
        except ImportError as e:
            print " `-- module not found.", e
        except KeyError as e:
            print " `-- tool not found.", e


class ProjectInfo( PyTool ):
    """
    """
    args = [
        _( "project_file", type="file", ext="json" )
    ]
    out = []

    def _init( self, *args, **kwargs ):
        self.project = read_project_file( self.project_file )

    def func( self ):
        if not self.project:
            return
        project_info( self.project )


class ProjectRun( PyTool ):
    """
    """
    args = [
        _( "project_file", type="file", ext="json" ),
        _( "parts", type="str", nargs="*", default=None,
            help="parts to consider, leave empty for all" ),
        _( "tools", type="str", nargs="*", default=None,
            help="tools to consider, leave empty for all" ),
        _( "print_args|pa", type="bool", default=False ),
        _( "simulate|sim", type="bool", default=False ),
        _( "verbose_tool|vt", type="bool", default=False ),
        _( "traceback|tr", type="bool", default=False ),

        # mdkit extension
        _( "analyze_only|ao", type="int", default=None ),
    ]
    out = []

    def _init( self, *args, **kwargs ):
        self.project = read_project_file( self.project_file )

    def func( self ):
        if not self.project:
            return

        def args_copy( args, defaults ):
            new_args = copy.deepcopy( args )
            new_args.update( defaults )
            if ( "__append__" in args and "__append__" in defaults ):
                new_args["__append__"] = copy.deepcopy(
                    args["__append__"]
                )
                new_args["__append__"].update(
                    defaults["__append__"]
                )
            return new_args

        def pre_process( tools ):
            tools_pp = collections.OrderedDict()
            for tid, t in tools.iteritems():
                if "__each__" in t:
                    each = t.pop("__each__")
                    each2 = None
                    each3 = None
                    if "__each__" in each:
                        each2 = each.pop("__each__")
                        if "__each__" in each2:
                            each3 = each2.pop("__each__")
                    for name, defaults in each.iteritems():
                        d = args_copy( t, defaults )
                        if each2:
                            for name2, defaults2 in each2.iteritems():
                                d2 = args_copy( d, defaults2 )
                                if each3:
                                    for name3, defaults3 in each3.iteritems():
                                        d3 = args_copy( d2, defaults3 )
                                        tools_pp[ tid + name + name2 + name3 ] = d3
                                else:
                                    tools_pp[ tid + name + name2 ] = d2
                        else:
                            tools_pp[ tid + name ] = d
                else:
                    tools_pp[ tid ] = t
            return copy.deepcopy( tools_pp )
        self.project["tools"] = pre_process( self.project["tools"] )

        parts_pp = collections.OrderedDict()
        for pid, p in self.project["parts"].iteritems():
            if "__parent__" in p:
                parent = p["__parent__"]
                del p["__parent__"]
                p2 = copy.deepcopy(
                    self.project["parts"][ parent ]
                )
                parts_pp[ pid ] = p2
            else:
                parts_pp[ pid ] = p
        self.project["parts"] = parts_pp

        for pid, p in self.project["parts"].iteritems():
            if p.get("__sub__"):
                self.project["parts"][pid]["__sub__"] = pre_process(
                    p["__sub__"]
                )

        for pid, p in self.project["parts"].iteritems():
            if p.get("__sub__"):
                sub_pp = collections.OrderedDict()
                for pid, sub in p["__sub__"].iteritems():
                    if "__range__" in sub or "__list__" in sub:
                        if "__range__" in sub:
                            r = sub.pop("__range__")
                            range_list = range( r[0], r[1] + 1 )
                        elif "__list__" in sub:
                            range_list = sub.pop("__list__")
                        for i in range_list:
                            sub_pp[ pid + str(i) ] = copy.deepcopy(
                                p["__sub__"][pid]
                            )
                            sub_pp[ pid + str(i) ][ "i" ] = i
                            for k, v in sub_pp[ pid + str(i) ].iteritems():
                                if isinstance( v, basestring ):
                                    v = v.format( i=i )
                                    sub_pp[ pid + str(i) ][k] = v
                        del p["__sub__"][pid]
                    else:
                        sub_pp[ pid ] = sub
                p["__sub__"] = sub_pp

        parts_pp = collections.OrderedDict()
        for pid, p in self.project["parts"].iteritems():
            parts_pp[pid] = p
            if p.get("__sub__"):
                sub = p.pop("__sub__")
                for sid, s in sub.iteritems():
                    parts_pp[ pid + "#" + sid ] = copy.deepcopy( p )
                    parts_pp[ pid + "#" + sid ].update( s )
                    parts_pp[ pid + "#" + sid ]["__dir__"] = os.path.join(
                        p["__dir__"], s["__dir__"]
                    )
        self.project["parts"] = parts_pp

        with open( "foo.json", "w" ) as fp:
            json.dump( self.project, fp, indent=4 )

        print "{:#^80}".format( " FILTER " )

        self.project["parts_all"] = copy.deepcopy( self.project["parts"] )
        for x, s in [ (self.parts, "parts"), (self.tools, "tools") ]:
            print "%s: %s" % (s.upper(), ", ".join(x) if x else "all")
            if x:
                self.project[s] = collections.OrderedDict([
                    (yid, y) for (yid, y) in self.project[s].iteritems()
                    if any([ re.match( "^" + xid + "$", yid ) for xid in x ])
                ])

        project_info( self.project )

        print "{:#^80}".format(
            " SIMULATE " if self.simulate else " RUN "
        )

        for tid, t in self.project["tools"].iteritems():
            print "ID: {id:<45}NAME: {__name__}".format( id=tid, **t )
            tool = get_tool( t )

            if "__summary__" in t:
                kwargs = copy.deepcopy( self.get_kwargs( tid, t, {} ) )

                print " `-- SUMMARY: %s" % ", ".join([ "all" ])

                parts = copy.deepcopy( self.project["parts"] )
                parts_all = copy.deepcopy( self.project["parts_all"] )
                if "__parts__" in kwargs:
                    parts = collections.OrderedDict([
                        (pid, parts_all[pid]) for pid in
                        kwargs["__parts__"] if pid in parts_all
                    ])

                if "__append__" in kwargs:
                    d = collections.defaultdict(list)
                    for pid, p in parts.iteritems():
                        for k, v in kwargs["__append__"].iteritems():
                            v = listify( v )
                            v2 = []
                            for _v in v:
                                if isinstance( _v, basestring ):
                                    _v = _v.format(
                                        dir=get_wd( self.project, p ),
                                        **self.get_kwargs( tid, t, p )
                                    )
                                v2.append( _v )
                            d[ k ] += v2
                    del kwargs["__append__"]
                    kwargs.update( d )

                for k, v in kwargs.iteritems():
                    if v == "$part_name_list":
                        kwargs[ k ] = parts.keys()
                    elif isinstance( v, basestring ):
                        kwargs[ k ] = v.format(
                            pid="|".join( parts.keys() ),
                            tid=tid, **kwargs
                        )

                wd = self.output_dir.format( tid=tid, **kwargs )

                if self.print_args:
                    print "      `-- WD: %s" % wd
                    print json.dumps( kwargs, indent=4 )

                if not self.simulate:
                    ret = self.call_tool( tool, kwargs, wd )
                    print "      `-- TOOL: %s" % ret

            else:
                for pid, p in self.project["parts"].iteritems():
                    print " `-- ID: {id:<20}DIR: {__dir__}".format(
                        id=pid, **p
                    )

                    kwargs = self.get_kwargs( tid, t, p )

                    for k, v in kwargs.iteritems():
                        if isinstance( v, basestring ):
                            if v.startswith("$") and "__variables__" in p:
                                for var, d in p["__variables__"].iteritems():
                                    if v == var:
                                        kwargs[ k ] = d
                                        break
                                else:
                                    print "%s not found in __variables__" % v
                            else:
                                kwargs[ k ] = v.format(
                                    pid=pid, tid=tid, **kwargs
                                )

                    wd = get_wd( self.project, p, pid=pid, tid=tid, **kwargs )

                    if "__variables__" in kwargs:
                        del kwargs["__variables__"]

                    if self.print_args:
                        print "      `-- WD: %s" % wd
                        print json.dumps( kwargs, indent=4 )

                    if not self.simulate:
                        ret = self.call_tool( tool, kwargs, wd )
                        print "      `-- TOOL: %s" % ret

    def get_kwargs( self, tid, t, p ):
        kwargs = {
            "analyze_only": self.analyze_only,
            "run": True,
            "verbose": self.verbose_tool,
        }
        kwargs.update( self.project.get("defaults", {}) )
        kwargs.update( t )
        kwargs.update( p )
        if "__tool_defaults__" in kwargs:
            kwargs.update( kwargs["__tool_defaults__"].get( tid, {} ) )
        copy.deepcopy( kwargs )
        if "__tool_defaults__" in kwargs:
            del kwargs["__tool_defaults__"]
        return kwargs

    def get_args( self, tool, kwargs ):
        args = []
        for a, params in tool.args.iteritems():
            if "default" not in params:
                args.append( kwargs.pop(a) )
        return args

    def call_tool( self, tool, kwargs, wd ):
        args = self.get_args( tool, kwargs )
        if not os.path.exists( wd ):
            os.makedirs( wd )
        with working_directory( wd ):
            try:
                ret = tool( *args, **kwargs )
            except Exception as e:
                if self.traceback:
                    import traceback
                    traceback.print_exc()
                ret = "Exception, %s" % e
        return ret
