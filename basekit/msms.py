#! /usr/bin/env python

from __future__ import with_statement
from __future__ import division

import numpy as np
import utils.numpdb as numpdb
import os
import string
import shutil
import collections

import utils
from utils import copy_dict, dir_walker, iter_stride
from utils.tool import _, _dir_init, CmdTool, ProviMixin
from utils.math import hclust
import json

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "msms" )
wfmesh_FILE = os.path.join( TMPL_DIR, "wfmesh.py" )
MSMS_CMD = "msms"
BABEL_CMD = "babel"


# newest open babel helps 2.2.3


"""
g Surface1
# Number of vertices: 30918
v -51.84846 50.86638 114.37491

# Number of normals: 21610
vn -.076 .994 -.082

# Number of faces: 64556
f 3//3 2//2 1//1

"""
def msms_to_obj(  ):
    pass



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
            # find errors
            if line.find( "ERROR Too many RS components" )!=-1:
                raise Exception( "too many RS components" )
            if line.find( "ERROR: find_first_rs_face" )!=-1:
                raise Exception( "find_first_rs_face" )
            # if line.find( "sphere_mange_arete: inconcistence" )!=-1:
            #     raise Exception( "sphere_mange_arete: inconcistence" )
            if ( line.find( "RS component" )!=-1 and 
                    line.find( "not found" )!=-1 ):
                raise Exception( "RS component not found" )
            if ( line.find( "input file" )!=-1 and 
                    line.find( "couldn't be open" )!=-1 ):
                raise Exception( "input file couldn't be open" )
            if line.find( 
                    "failed after 5 attemps to compute the surface" )!=-1:
                raise Exception( 
                    "failed after 5 attemps to compute the surface" )
            
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
            if line.find( "ERROR: is_point_in_contact_face" )!=-1:
                pass
                # raise Exception( "is_point_in_contact_face" )
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



def parse_msms_area( area_file, max_atomno=None ):
    with open( area_file, "r" ) as fp:
        header = fp.next().split()[1:]
        ses_list = [ [] for i in range(len(header)) ]
        sas_list = [ [] for i in range(len(header)) ]
        for line in fp:
            if not line.strip():
                continue
            ls = line.split()
            atomno = int(ls[0]) + 1
            if max_atomno and atomno>max_atomno:
                break
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


def make_nrhole_pdb( pdb_input, holes, nh_file, mean_file, pymol_file, obj_list):
    npdb = numpdb.NumPdb( pdb_input )
    try:
        sele2={'record':'HETATM'}
        last_hetatm=int(npdb.get('atomno', **sele2)[-1])
    except:
        sele2={'record':'ATOM  '}
        last_hetatm=int(npdb.get('atomno', **sele2)[-1])
    fi=open(nh_file, 'w')
    preline=""
    mean_dct={}
    neighbours={}
    with open(pdb_input, 'r') as fp:
        for line in fp:
            if (preline!=line[0:6] and preline=="HETATM")  or line[0:3]=="TER" or line[0:3]=="END":
                for hno, holeneighbours in enumerate(holes["ses"]):
                    if holeneighbours!=[] and hno!=0:
                        chain=''
                        xyz_list=[]
                        atom_list=[]
                        for hn in holeneighbours:
                            try:
                                atomno=hn
                                sele={'atomno':atomno}
                                resno=npdb.get('resno', **sele)[0]
                                sele={'atomno':atomno}
                                chain=npdb.get('chain', **sele)[0]
                                sele={'atomno':atomno, 'chain':chain }
                                xyz=npdb.get('xyz', **sele)[0]
                                xyz_list.append(xyz)
                                atom_list.append([atomno, int(resno), chain])
                            except:
                                break
                        neighbours[hno]=atom_list
                        xyz2=xyz_list
                        coords = np.mean(xyz_list, axis=0)
                        xyz_list_dist=[]
                        for co in xyz2:
                            xyz_list_dist.append(np.sqrt(np.sum((coords-co)**2)))
                        mean=np.mean(xyz_list_dist, axis=0)
                        std=np.std(xyz_list_dist)
                        natom={
                            "record":"HETATM","atomno": last_hetatm+hno, "atomname": " CA ",
                            "resname": "NEH", "chain": chain, "resno": hno,
                            "x": coords[0], "y": coords[1], "z": coords[2], "bfac": std, "element": " C"
                        }
                        mean_dct[last_hetatm+hno]=mean
                        fi.write(numpdb.pdb_line( natom ))
                fi.write(line)
            else:
                fi.write(line)
            preline=line[0:6]
    with open( mean_file, "w" ) as fp:
        json.dump( mean_dct, fp )
    return json.dumps(mean_dct), json.dumps(neighbours), obj_list
    #make_nrholes_pymol(std_file, mean_dct, pymol_file, neighbours, nh_file, obj_list)

def f_v_to_obj(face_file, vert_file, mesh_file, flip=False):
    def writer(f, fp, pre):
        tmp=True
        with open(f, 'r') as fp2:
            for line in fp2:
                if line.startswith('#'):
                    fp.write(line)
                elif tmp:
                    tmp=False
                else:
                    splitted=line.split()
                    if pre=='f':
                        fp.write(pre+' '+splitted[0]+' '+splitted[2]+' '+splitted[1]+'\n') 
                    else:
                        fp.write(pre+' '+splitted[0]+' '+splitted[1]+' '+splitted[2]+'\n')    
    with open(mesh_file, 'w') as fp:
        writer(vert_file, fp, 'v')
        writer(face_file, fp, 'f')

class Msms( CmdTool, ProviMixin ):
    """A wrapper around the MSMS program."""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "probe_radius", type="float", range=[0.1, 5], step=0.1, 
            default=1.5 ),
        _( "density", type="float", range=[0.5, 10], step=0.5, default=1.0 ),
        _( "hdensity", type="float", range=[1.0, 20], step=1.0, default=3.0 ),
        _( "all_components", type="bool", default=False ),
        _( "no_area", type="bool", default=False ),
        _( "envelope", type="float", range=[0.0, 10], step=0.1, 
            default=0 ),
        _( "envelope_hclust", type="str", default="", 
            options=[ "", "ward", "average" ], help="'', average, ward" ),
        _( "atom_radius_add", type="float", default=0 ),
        _( "getinfo", type="bool", default=False ),
        _( "get_nrholes|gh", type="bool", default=False ),
    ]
    out = [
        _( "area_file", file="area.area" ),
        _( "face_file", file="tri_surface.face" ),
        _( "vert_file", file="tri_surface.vert" ),
        _( "surf_file", file="surf.xyzr" ),
        _( "surf_file2", file="surf2.xyz" ),
        _( "surf_file3", file="surf3.xyz" ),
        _( "surf_file4", file="surf4.xyz", optional=True ),
        _( "provi_file", file="msms.provi" ),
        _( "info_file", file="info_file.txt", optional=True ),
        _( "nh_file", file="{pdb_file.stem}_nh.pdb", optional=True ),
        _( "pymol_file", file="pymol_settings.py", optional=True ),
        _( "mean_file", file="{pdb_file.stem}_mean.json", optional=True ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "msms.provi"
    pymol_tmpl = "pymol_settings.py"
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
            if self.atom_radius_add:
                with open( self.pdb2xyzr.xyzr_file, "r+" ) as fp:
                    lines = fp.readlines()
                    fp.seek(0)
                    fp.truncate()
                    for l in lines:
                        x, y, z, r = l.split()
                        l = "%s\t%s\t%s\t%0.2f\n" % (
                            x, y, z, float(r) + self.atom_radius_add
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
        if self.getinfo:
            self.get_infos()
        v = "tri_surface(_?[0-9]*)\.(vert)"
        obj_list=[]
        for m, vertfile in dir_walker( self.output_dir, v ):
            facefile = vertfile[:-4]+'face'
            meshfile = vertfile[:-4]+'obj'
            num=meshfile.split('_')[-1].split('.')[0]
            try:
                obj_list.append([int(num), meshfile])
                f_v_to_obj(facefile, vertfile, meshfile, flip=True)
            except ValueError:
                f_v_to_obj(facefile, vertfile, meshfile, flip=False)
# get the nrholes and the pymol script
        if self.get_nrholes:
            area = parse_msms_area( self.area_file)
            mean_dct, neighbours, obj_list = make_nrhole_pdb( self.pdb_file, area, self.nh_file, self.mean_file, self.pymol_file, obj_list )
            values_dict={'neighbours':neighbours, 'obj_list':obj_list, 'mean_dct':mean_dct, 'nh_file':self.nh_file, 'TMPL_DIR':self.tmpl_dir}
            self._make_file_from_tmpl(self.pymol_tmpl, **values_dict)
    def get_infos( self ):
        def make_list(area):
            sas_cav={}
            for cav_no, cav in enumerate(area):
                for atom in cav:
                    if sas_cav.has_key(cav_no):
                        try:
                            sas_cav[cav_no].append(res_list[atom])
                        except:
                            pass
                    else:
                        try:
                            sas_cav[cav_no] = [res_list[atom]]
                        except:
                            pass
            return sas_cav
        comps = self.get_components()
        area = parse_msms_area( self.area_file)
        res_list={}
        with open(self.pdb_file, 'r') as fp:
            for line in fp:
                if line.startswith("ATOM"):
                    res_list[int(line[6:11])] = (line[21:22], int(line[22:26]))
        ses_cav=make_list(area["ses"])
        result_list=[["Cav", "ses_area", "ses_vol", "ses_neighbour_residue_list"]]
        for cav in range(len(comps)):
            r=sorted(list(set(ses_cav[cav])))
            result_list.append([comps[cav][0], comps[cav][1], comps[cav][2], r])
        with open( self.info_file, "w" ) as fp:
            json.dump( result_list, fp, indent=4 )
    def components_provi( self, color="", translucent=0.0, relpath=None,
                            max_atomno=None ):
        if self.all_components:
            comps = self.get_components()
            if max_atomno:
                area = self.get_area( max_atomno=max_atomno )
                comps = [ c for c in comps if area["ses"][ c[0] ] ]
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
    def get_area( self, max_atomno=None ):
        return parse_msms_area( self.area_file, max_atomno=max_atomno )
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
        else:
            xyz_file = utils.path.mod( self.xyzr_file, ext="xyz" )
            with open( self.xyzr_file, "r" ) as fp:
                lines = fp.readlines()
                with open( xyz_file, "w" ) as fp_out:
                    fp_out.write( "%i\ncomment\n" % len( lines ) )
                    for d in lines:
                        ds = d.strip().split()
                        if not len(ds):
                            continue
                        l = "F\t%0.3f\t%0.3f\t%0.3f\n" % (
                            float(ds[0]), float(ds[1]), float(ds[2])
                        )
                        fp_out.write( l )
