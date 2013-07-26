#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division


import re
import sys
import os
import shutil
import collections

import utils.path
from utils import copy_dict, dir_walker, iter_stride
from utils.tool import _, _dir_init, CmdTool, ProviMixin

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "msms" )

MSMS_CMD = "msms"
BABEL_CMD = "babel"


MsmsComponent = collections.namedtuple( "MsmsComponent", [
    "id", "SES_volume", "SES_area"
])
def parse_msms_log( msms_log ):
    components = [] #collections.OrderedDict()
    with open( msms_log, "r" ) as fp:
        for line in fp:
            if ( line.startswith("ANALYTICAL SURFACE AREA") or
                    line.startswith("TRIANGULATION") or
                    line.startswith("NUMERICAL VOLUMES AND AREA") ):
                break
        fp.next()
        for line in fp:
            if ( line.startswith("TRIANGULATION") or
                    line.startswith("NUMERICAL VOLUMES AND AREA") ):
                break
            else:
                pass
                # print line.split()
        for line in fp:
            if line.startswith("NUMERICAL VOLUMES AND AREA"):
                break
        fp.next()
        for line in fp:
            if line.strip().startswith("Total"):
                break
            else:
                ls = line.split()
                components.append(( 
                    int(ls[0]),     # id
                    float(ls[2]),   # SES volume
                    float(ls[3])    # SES area
                ))
    return components



def parse_msms_area( area_file ):
    with open( area_file, "r" ) as fp:
        header = fp.next().split()[1:]
        ses_list = [ [] for i in range(len(header)) ]
        sas_list = [ [] for i in range(len(header)) ]
        for line in fp:
            if not line.strip():
                continue
            ls = line.split()
            atomno = int(ls[0]) + 1
            atom_area = map( float, ls[1:] )
            for i, a in enumerate( iter_stride( atom_area ) ):
                ses, sas = a
                if ses>0:
                    ses_list[i].append( atomno )
                if sas>0:
                    sas_list[i].append( atomno )
    return {
        "ses": ses_list,
        "sas": sas_list
    }




class Msms( CmdTool, ProviMixin ):
    """A wrapper around the MSMS program."""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "probe_radius", type="float", range=[0.1, 5], step=0.1, 
            default=1.5 ),
        _( "density", type="float", range=[0.5, 10], step=0.5, default=1.0 ),
        _( "hdensity", type="float", range=[1.0, 20], step=1.0, default=3.0 ),
        _( "all_components", type="checkbox", default=False ),
        _( "no_area", type="checkbox", default=False ),
    ]
    out = [
        _( "area_file", file="area.area" ),
        _( "face_file", file="tri_surface.face" ),
        _( "vert_file", file="tri_surface.vert" ),
        _( "provi_file", file="msms.provi" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "msms.provi"
    def _init( self, *args, **kwargs ):
        self.pdb2xyzr = Pdb2xyzr( 
            self.pdb_file, **copy_dict( kwargs, run=False ) 
        )
        self.cmd = [ 
            MSMS_CMD, "-if", self.pdb2xyzr.xyzr_file, 
            "-probe_radius", self.probe_radius,
            "-af", "area", "-of", "tri_surface", 
            "-density", self.density,
            "-hdensity", self.hdensity
        ]
        if self.all_components:
            self.cmd.append( "-all_components" )
        if self.no_area:
            self.cmd.append( "-no_area" )
        self.output_files.extend( self.pdb2xyzr.output_files )
    def _pre_exec( self ):
        self.pdb2xyzr()
    def _post_exec( self ):
        if self.all_components:
            p = "tri_surface_([0-9]+)\.vert"
            components = []
            for m, filepath in dir_walker( self.output_dir, p ):
                components.append(
                    '\t{\n\t\t"filename": "%s"\n\t}' % self.relpath( filepath )
                )
            components_data = ",\n" + ",\n".join( components )
        else:
            components_data = ""

        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            vert_file=self.relpath( self.vert_file ),
            components_data=components_data
        )
    def get_components( self ):
        return parse_msms_log( self.stdout_file )
    def get_area( self ):
        return parse_msms_area( self.area_file )




class Pdb2xyzr( CmdTool ):
    """A pdb to xyzr format converter based on OpenBabel."""
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    out = [
        _( "pdb_prep_file", file="{pdb_file.stem}_prep.pdb" ),
        _( "xyzr_file", file="{pdb_file.stem}.xyzr" )
    ]
    def _init( self, *args, **kwargs ):
        self.cmd = [ 
            BABEL_CMD, '-i', 'pdb', self.pdb_prep_file,
            '-o', 'msms', self.xyzr_file 
        ]
        self.output_files = [ self.pdb_prep_file, self.xyzr_file ]
    def _pre_exec( self ):
        shutil.copy( self.pdb_file, self.pdb_prep_file )


