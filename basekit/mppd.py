from __future__ import with_statement
from __future__ import division

import math
import collections

import numpy as np

import utils.path
from utils import copy_dict
from utils.math import point_plane_dist
from utils.tool import _, _dir_init, PyTool, ProviMixin
from utils.numpdb import NumPdb
from opm import Opm
from dowser import DowserRepeat
from voronoia import Voronoia, HOLE_NOT_FILLED, HOLE_PARTLY_FILLED
from hbexplore import HBexplore
from msms import Msms


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "mppd" )



class MppdPipeline( PyTool, ProviMixin ):
    """The MPPD pipeline"""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "analyze_only|ao", type="checkbox", default=False ),
        _( "variants", type="str", nargs="*", default=[] ),
        _( "tools", type="str", nargs="*", default=[] )
    ]
    out = [
        _( "processed_pdb", file="{pdb_file.stem}_proc.pdb" ),
        _( "no_water_file", file="{pdb_file.stem}_nowat.pdb" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "mppd.provi"
    def _init( self, *args, **kwargs ):
        self.pdb_id = utils.path.stem( self.pdb_file )
        self.probe_radius = 1.4
        self.opm = Opm(
            self.pdb_id,
            **copy_dict( kwargs, run=False, 
                output_dir=self.subdir("opm") )
        )
        self.dowser = DowserRepeat(
            self.processed_pdb,
            **copy_dict( kwargs, run=False, 
                output_dir=self.subdir("dowser") )
        )
        self.water_variants = [
            ( "org", self.processed_pdb ),
            ( "non", self.no_water_file ),
            ( "dow", self.dowser.dowser_file )
        ]
        self.tool_list = [
            ( "voronoia", Voronoia, {} ),
            ( "hbexplore", HBexplore, {} ),
            ( "msms", Msms, { 
                "all_components": True, 
                "probe_radius": self.probe_radius 
            })
        ]
        for suffix, pdb_file in self.water_variants:
            for prefix, tool, tool_kwargs in self.tool_list:
                name = "%s_%s" % ( prefix, suffix )
                self.__dict__[ name ] = tool(
                    pdb_file, **copy_dict( kwargs, run=False, 
                        output_dir=self.subdir( name ), **tool_kwargs )
                )
                self.output_files += self.__dict__[ name ].output_files
    def _pre_exec( self ):
        self.make_processed_pdb()
        self.make_nowat_pdb()
    def func( self ):
        if not self.analyze_only:
            water_variants = [ 
                v for v in self.water_variants 
                if v[0] in self.variants or not self.variants 
            ]
            tool_list = [ 
                t for t in self.tool_list 
                if t[0] in self.tools or not self.tools 
            ]
            for name in [ "opm", "dowser" ]:
                if name in self.tools or not self.tools:
                    self.__dict__[ name ]()
            for suffix, pdb_file in water_variants:
                for prefix, tool, tool_kwargs in tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ]()
        self.stats()
    def _post_exec( self ):
        self._make_provi_file(
            family="?",
            probe_radius=self.probe_radius,
            pdb_id=utils.path.stem( self.pdb_file ),
            pdb_file=self.relpath( self.processed_pdb ),
            dowser_file=self.relpath( self.dowser.dowser_file ),
            vol_file_org=self.relpath( self.voronoia_org.vol_file ),
            hbx_file_org=self.relpath( self.hbexplore_org.hbx_file ),
            vol_file_dow=self.relpath( self.voronoia_dow.vol_file ),
            hbx_file_dow=self.relpath( self.hbexplore_dow.hbx_file )
        )
    def make_processed_pdb( self ):
        npdb = NumPdb( self.opm.processed_file )
        npdb.write2( self.processed_pdb )
    def make_nowat_pdb( self ):
        with open( self.no_water_file, "w" ) as fp:
            with open( self.processed_pdb, "r" ) as fp_pdb:
                for line in fp_pdb:
                    if ( line[0:6] in ["ATOM  ", "HETATM"] and
                            line[17:20]!="HOH" ):
                        fp.write( line )
    def get_npdb_dicts( self ):
        mplanes = self.opm.get_planes()
        dist = abs( point_plane_dist( mplanes[1][0], mplanes[0] ) )
        npdb_dict = {}
        npdb_tm_dict = {}
        for suffix, pdb_file in self.water_variants:
            npdb = NumPdb( pdb_file )
            npdb_dict[ suffix ] = npdb
            sele = npdb.sele()
            i = 0
            # create transmembrane region selection
            for numa in npdb.iter_resno( incomplete=True ):
                flag = True
                for c in numa._coords:
                    d1 = abs( point_plane_dist( c, mplanes[0] ) )
                    d2 = abs( point_plane_dist( c, mplanes[1] ) )
                    if d1<dist and d2<dist:
                        flag = False
                        break
                for a in numa._atoms:
                    sele[i] = flag
                    i += 1
            np.logical_not( sele, sele )
            npdb_tm_dict[ suffix ] = npdb.copy( sele=sele )
            npdb.copy( sele=sele ).write( self.outpath( "tm_region.pdb" ) )
        return npdb_dict, npdb_tm_dict
    def stats( self ):
        npdb_dict, npdb_tm_dict = self.get_npdb_dicts()
        mplanes = self.opm.get_planes()
        dist = abs( point_plane_dist( mplanes[1][0], mplanes[0] ) )
        segments = [ ("ALL", npdb_dict), ("TM", npdb_tm_dict) ]
        pr = self.probe_radius
        mppd_records = []
        for suffix, pdb_file in self.water_variants:
            hbx = self.__dict__[ "hbexplore_%s" % suffix ]
            voro = self.__dict__[ "voronoia_%s" % suffix ]
            msms = self.__dict__[ "msms_%s" % suffix ]
            for seg, npdb_d in segments:
                mppd_records.append( MppdRecord(*(
                    ( self.pdb_id, suffix, seg, dist ) +
                    count_water( npdb_d[ suffix ] ) +
                    ( None, ) +
                    count_hbonds( npdb_dict[ suffix ], hbx ) +
                    count_holes_voro( npdb_dict[ suffix ], voro ) +
                    count_holes_msms( npdb_tm_dict[ suffix ], msms, pr )
                )))
        for r in mppd_records:
            r.info()




_MppdRecord = collections.namedtuple( "_MppdRecord", [
    "pdb_id", "source", "segment", "mplanes_dist", 
    "water_count", "residue_count", "hetero_count", 
    "chain_count", "identical_chains",
    "hbx_protein_water", "hbx_water_water",
    "voro_not_filled", "voro_partly_filled",
    "msms", "msms_ses", "msms_gt_water", "msms_gt_water_ses"
])
class MppdRecord( _MppdRecord ):
    def info( self ):
        print "### %s [ %s, %s ]" % ( 
            self.pdb_id, self.source, self.segment 
        )
        print "distance: %0.2f" % self.mplanes_dist
        print "counts: %i (water), %i (aa), %i (hetero)" % (
            self.water_count, self.residue_count, self.hetero_count
        )
        print "chains: %i, identical: %s" % (
            self.chain_count, "??? (TODO)"
        )
        print "hbx: %i (protein-water), %i (water-water)" % (
            self.hbx_protein_water, self.hbx_water_water
        )
        print "voro: %i (not filled), %i (partly filled)" % (
            self.voro_not_filled, self.voro_partly_filled
        )
        print "msms: %i, %0.2f A^3; %i (SES>H2O), %0.2f A^3" % (
            self.msms, self.msms_ses, 
            self.msms_gt_water, self.msms_gt_water_ses
        )
        print ""




def count_water( npdb ):
    water_count = 0
    aa_count = 0
    hetero_count = 0
    chain_count = len( np.unique( 
        npdb.copy(resname="HOH", invert=True)["chain"] 
    ))
    for numa in npdb._iter_resno():
        if numa[0]["resname"]=="HOH" and numa[0]["atomname"]==" O  ":
            water_count += 1
        elif numa[0]["record"]=="ATOM  ":
            aa_count += 1
        elif numa[0]["record"]=="HETATM":
            hetero_count += 1
    return ( water_count, aa_count, hetero_count, chain_count )


def count_hbonds( npdb, hbexplore ):
    hb_count = 0
    hb_count2 = 0
    for hb in hbexplore.get_hbonds():
        # check if both atoms are in the npdb
        na1 = npdb.copy( resno=hb.resno1, chain=hb.chain1 )
        na2 = npdb.copy( resno=hb.resno2, chain=hb.chain2 )
        if len(na1) and len(na2):
            if hb.resname1=="HOH" and hb.resname2=="HOH":
                hb_count2 += 1
            elif hb.resname1=="HOH" or hb.resname2=="HOH":
                hb_count += 1
    return ( hb_count, hb_count2 )


def count_holes_voro( npdb, voronoia ):
    vol = voronoia.get_vol()
    not_filled = 0
    partly_filled = 0
    for hole in vol[ "holes" ]:
        for nb in hole.neighbours:
            if not len( npdb.copy( atomno=nb.atomno ) ):
                break
        else:
            if hole.type==HOLE_NOT_FILLED:
                not_filled += 1
            elif hole.type==HOLE_PARTLY_FILLED:
                partly_filled += 1
    return ( not_filled, partly_filled )


def count_holes_msms( npdb, msms, probe_radius ):
    area = msms.get_area()
    components = msms.get_components()
    components0 = []
    for c in components:
        for nb_atomno in area["ses"][ c[0] ]:
            if not len( npdb.copy( atomno=nb_atomno ) ):
                break
        else:
            components0.append( c )
    components2 = [
        c for c in components0 if c[1] < 0
    ]
    water_vol = (4.0/3.0) * math.pi * (probe_radius**3)
    components3 = [
        c for c in components0 if c[1] < -water_vol
    ]
    return (
        len( components2 ), 
        sum([ c[1] for c in components2]),
        len( components3 ), 
        sum([ c[1] for c in components3])
    )


