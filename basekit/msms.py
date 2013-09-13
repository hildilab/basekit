#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division

import numpy as np
import scipy.cluster
import scipy.spatial

import os
import string
import shutil
import collections

import utils
from utils import copy_dict, dir_walker, iter_stride
from utils.tool import _, _dir_init, CmdTool, ProviMixin
from utils.math import hclust

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
            if line.find( "ERROR Too many RS components" )!=-1:
                raise Exception( "too many RS components" )
            if line.find( "ERROR: find_first_rs_face" )!=-1:
                raise Exception( "find_first_rs_face" )
            if line.find( "sphere_mange_arete: inconcistence" )!=-1:
                raise Exception( "sphere_mange_arete: inconcistence" )
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
            for i, a in enumerate( iter_stride( atom_area, n=2 ) ):
                ses, sas = a
                if ses>0:
                    ses_list[i].append( atomno )
                if sas>0:
                    sas_list[i].append( atomno )
    return {
        "ses": ses_list,
        "sas": sas_list
    }



def parse_msms_vert( vert_file ):
    vert_list = []
    with open( vert_file, "r" ) as fp:
        header = fp.next() + fp.next() + fp.next()
        for line in fp:
            ls = line.split()
            vert_list.append((
                float( ls[0] ), float( ls[1] ), float( ls[2] ),
                float( ls[3] ), float( ls[4] ), float( ls[5] ),
                int( ls[6] ), int( ls[7] ), int( ls[8] )
            ))
    types = [
        ('x', np.float), ('y', np.float), ('z', np.float),
        ('nx', np.float), ('ny', np.float), ('nz', np.float),
        ('face_no', np.int),
        ('closest_sphere', np.int),
        ('face_type', np.int)
    ]
    return np.array( vert_list, dtype=types )



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
        _( "envelope", type="float", range=[0.0, 10], step=0.1, 
            default=0 ),
        _( "envelope_hclust", type="str", default="average" 
            options=[ "", "ward", "average" ], help="'', average, ward" ),
    ]
    out = [
        _( "area_file", file="area.area" ),
        _( "face_file", file="tri_surface.face" ),
        _( "vert_file", file="tri_surface.vert" ),
        _( "surf_file", file="surf.xyzr" ),
        _( "surf_file2", file="surf2.xyz" ),
        _( "surf_file3", file="surf3.xyz" ),
        _( "surf_file4", file="surf4.xyz", optional=True ),
        _( "provi_file", file="msms.provi" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "msms.provi"
    def _init( self, *args, **kwargs ):
        self.pdb2xyzr = Pdb2xyzr( 
            self.pdb_file, **copy_dict( kwargs, run=False, verbose=False ) 
        )
        self.output_files = self.pdb2xyzr.output_files + self.output_files
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
        if self.envelope:
            self.envelope_pdb = self.outpath( "envelope.pdb" )
            self.envelope_msms = Msms(
                self.pdb_file, **copy_dict( kwargs, run=False,
                    envelope=0.0, density=0.3, hdensity=0.3,
                    probe_radius=self.envelope, all_components=False,
                    output_dir=self.subdir( "envelope" ) )
            )
            self.output_files += self.envelope_msms.output_files
    def _pre_exec( self ):
        p = "tri_surface_([0-9]+)\.(vert|face)"
        for m, filepath in dir_walker( self.output_dir, p ):
            utils.path.remove( filepath )
        self.pdb2xyzr()
        if self.envelope:
            # calc envelope; append to self.pdb2xyzr.xyzr_file 
            self.envelope_msms()
            vert = self.envelope_msms.get_vert( filt=True )
            pr = self.envelope_msms.probe_radius
            if self.envelope_hclust:
                coords = np.array([ [ d['x'], d['y'], d['z'] ] for d in vert ])
                clust = hclust(
                    coords, 2*pr, method=self.envelope_hclust
                )
                print len( coords ), len( clust )
                with open( self.pdb2xyzr.xyzr_file, "a" ) as fp:
                    for d in clust.itervalues():
                        d3 = np.array( d ).mean( axis=0 )
                        l = "%0.3f\t%0.3f\t%0.3f\t%0.2f\n" % (
                            d3[0], d3[1], d3[2], pr
                        )
                        fp.write( l )
            else:
                with open( self.pdb2xyzr.xyzr_file, "a" ) as fp:
                    dct = {}
                    for d in vert:
                        coords = ( d['x'], d['y'], d['z'] )
                        if coords in dct:
                            continue
                        dct[ coords ] = True
                        l = "%0.3f\t%0.3f\t%0.3f\t%0.2f\n" % (
                            d['x'], d['y'], d['z'], pr
                        )
                        fp.write( l )
    def _post_exec( self ):
        self.get_components()
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            vert_file=self.relpath( self.vert_file ),
            color="",
            components_data=self.components_provi()
        )
        vert = self.get_vert()
        with open( self.surf_file, "w" ) as fp:
            for d in vert:
                l = "%0.3f\t%0.3f\t%0.3f\t%0.2f\n" % (
                    d['x'], d['y'], d['z'], self.probe_radius
                )
                fp.write( l )
        vert2 = filter( lambda x: x[-3]<0, vert )
        with open( self.surf_file2, "w" ) as fp:
            fp.write( "%i\ncomment\n" % len( vert2 ) )
            for d in vert2:
                l = "F\t%0.3f\t%0.3f\t%0.3f\n" % (
                    d['x'], d['y'], d['z']
                )
                fp.write( l )
        with open( self.surf_file3, "w" ) as fp:
            fp.write( "%i\ncomment\n" % len( vert ) )
            for d in vert:
                l = "F\t%0.3f\t%0.3f\t%0.3f\n" % (
                    d['x'], d['y'], d['z']
                )
                fp.write( l )
        if self.envelope_hclust:
            coords = np.array([ [ d['x'], d['y'], d['z'] ] for d in vert2 ])
            clust = hclust(
                coords, 2*self.probe_radius, method=self.envelope_hclust
            )
            with open( self.surf_file4, "w" ) as fp:
                fp.write( "%i\ncomment\n" % len( clust ) )
                for d in clust.itervalues():
                    d3 = np.array( d ).mean( axis=0 )
                    l = "F\t%0.3f\t%0.3f\t%0.3f\n" % (
                        d3[0], d3[1], d3[2]
                    )
                    fp.write( l )
    def components_provi( self, color="", translucent=0.0, relpath=None ):
        if self.all_components:
            comps = self.get_components()
            if not relpath:
                relpath = self.relpath
            with open( self.datapath( "_component.json" ), "r" ) as fp:
                components_tmpl = fp.read()
            p = "tri_surface_([0-9]+)\.vert"
            components = []
            for m, filepath in dir_walker( self.output_dir, p ):
                cno = int( m.group(1) )
                if cno >= len(comps) or comps[ cno ][1] > 0:
                    continue
                components.append( 
                    string.Template( components_tmpl ).substitute( 
                        vert_file=relpath( filepath ),
                        insideout="true",
                        color=color,
                        translucent=translucent
                    )
                )
            if components:
                return ",\n" + ",\n".join( components )
        return ""
    def component_files( self ):
        files = []
        p = "tri_surface_([0-9]+)\.(vert|face)"
        for m, filepath in dir_walker( self.output_dir, p ):
            files.append( os.path.abspath( filepath ) )
        return files
    def get_components( self ):
        return parse_msms_log( self.stdout_file )
    def get_area( self ):
        return parse_msms_area( self.area_file )
    def get_vert( self, filt=False ):
        vert = parse_msms_vert( self.vert_file )
        if filt:
            return filter( lambda x: x[-3]<0, vert )
        else:
            return vert




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
    def _post_exec( self ):
        if os.path.getsize( self.xyzr_file )==0:
            utils.path.remove( self.xyzr_file )

