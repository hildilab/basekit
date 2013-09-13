from __future__ import with_statement
from __future__ import division

import os
import re
import math
import json
import shutil
import collections

import numpy as np
import scipy.spatial
from scipy.stats.stats import pearsonr

# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

from utils import copy_dict, try_float, try_div
from utils.math import point_plane_dist
from utils.tool import (
    _, _dir_init, PyTool, RecordsMixin, ParallelMixin, ProviMixin,
    JsonBackend, SqliteBackend
)
from utils.numpdb import NumPdb, NumAtoms
from opm import Opm, OpmInfo
from dowser import DowserRepeat
from voronoia import Voronoia, HOLE_NOT_FILLED, HOLE_PARTLY_FILLED
from hbexplore import HBexplore
from msms import Msms
from mpstruc import MpstrucInfo
from pdb import PdbInfo


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "mppd" )



_MppdRecord = collections.namedtuple( "_MppdRecord", [
    "pdb_id", "source", "segment", "mplanes_dist", 
    "water_count", "residue_count", "hetero_count", 
    "chain_count", "identical_chains",
    "hbx_protein_water", "hbx_water_water", "hbx_water_ligand",
    "voro_not_filled", "voro_partly_filled", 
    "packdens_all", "packdens_protein", "packdens_water", "packdens_hetero",
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
        print "hbx: %i (protein-water), %i (water-water), %i (water-ligand)" % (
            self.hbx_protein_water, self.hbx_water_water, self.hbx_water_ligand
        )
        print "voro: %i (not filled), %i (partly filled)" % (
            self.voro_not_filled, self.voro_partly_filled
        )
        print "pd: %0.3f (all), %0.3f (aa), %0.3f (water), %0.3f (hetero)" % (
            self.packdens_all, self.packdens_protein, 
            self.packdens_water, self.packdens_hetero
        )
        print "msms: %i, %0.2f A^3; %i (SES>H2O), %0.2f A^3" % (
            self.msms, self.msms_ses, 
            self.msms_gt_water, self.msms_gt_water_ses
        )
        print ""


MppdDbRecord = collections.namedtuple( "MppdDbRecord", [
    "PDBID", "HEADER", "RESOLUTION", "EXPDTA", "KEYWDS",
    "OPMSpecies", "OPMFamily", "OPMSuperfamily",
    "swname", "swfamily", 
    "topic"
])



with open( os.path.join( TMPL_DIR, "pdb_dic_swfam.json" ), "r" ) as fp:
    PDB_DIC = json.load( fp )



def get_tree( coords ):
    if len( coords )==0:
        coords = np.array([[ np.inf, np.inf, np.inf ]])
    return scipy.spatial.KDTree( coords )


class MppdPipeline( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """The MPPD pipeline"""
    args = [
        _( "pdb_input", type="str" ),
        _( "probe_radius|pr", type="float", range=[0.1, 5],
            step=0.1, default=1.4, help="coulombic probe radius" ),
        _( "vdw_probe_radius|vpr", type="float", range=[0.1, 5], 
            step=0.1, default=1.7 ),
        _( "analyze_only|ao", type="checkbox", default=False ),
        _( "check_only|co", type="checkbox", default=False ),
        _( "variants", type="str", nargs="*", default=[], 
            help="a '!' as the first arg negates the list" ),
        _( "tools", type="str", nargs="*", default=[],
            help="a '!' as the first arg negates the list" ),
        _( "extract", type="str", default="" ),
        _( "figures|fig", type="checkbox", default=False ),
        _( "database|db", type="checkbox", default=False ),
        # msms tweaks
        _( "envelope_hclust", type="str", default="average",
            options=[ "", "ward", "average" ], help="'', average, ward" ),
    ]
    out = [
        _( "processed_pdb", file="proc.pdb" ),
        _( "no_water_file", file="nowat.pdb" ),
        _( "dowser_dry_pdb", file="dowser_dry.pdb" ),
        _( "original_dry_pdb", file="original_dry.pdb" ),
        _( "final_pdb", file="final.pdb" ),
        _( "stats_file", file="stats.json" ),
        _( "stats2_file", file="stats2.json" ),
        _( "info_file", file="info.json" ),
    ]
    RecordsClass = MppdRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "mppd.provi"
    pdb_dic = PDB_DIC
    def _init( self, *args, **kwargs ):
        self._init_records( None, **kwargs )
        self._init_parallel( self.pdb_input, **kwargs )

        if not self.parallel:
            self.pdb_id = self.pdb_input
            self.npdb_features = {
                "sstruc": False,
                "phi_psi": False
            }

            self.opm = Opm(
                self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("opm") )
            )
            self.output_files = self.opm.output_files
            self.output_files += [ self.processed_pdb, self.no_water_file ]
            self.dowser = DowserRepeat(
                self.processed_pdb,
                **copy_dict( kwargs, run=False, alt='x',
                    output_dir=self.subdir("dowser") )
            )
            self.output_files += self.dowser.output_files
            msms_kwargs = { 
                "all_components": True,
                "density": 1.0, "hdensity": 3.0,
                "envelope": self.probe_radius * 2,
                "envelope_hclust": self.msms_envelope_hclust
            }
            self.msms0 = Msms(
                self.no_water_file, **copy_dict( kwargs, run=False, 
                output_dir=self.subdir( "msms0" ), 
                **copy_dict( msms_kwargs, probe_radius=self.vdw_probe_radius)
            ))
            self.output_files += self.msms0.output_files
            self.output_files += [ self.original_dry_pdb, self.final_pdb ]

            self.pdb_info = PdbInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("pdb_info") ) )
            self.opm_info = OpmInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("opm_info") ) )
            self.mpstruc_info = MpstrucInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("mpstruc_info") ) )
            
            self.output_files += [ self.info_file ]

            self.water_variants = [
                ( "non", self.no_water_file ),
                ( "org", self.original_dry_pdb ),
                ( "dow", self.dowser_dry_pdb ),
                ( "fin", self.final_pdb ),
            ]
            
            self.tool_list = [
                ( "voronoia", Voronoia, { "ex":0.2 } ),
                ( "hbexplore", HBexplore, {} ),
                ( "msms_vdw", Msms, copy_dict( msms_kwargs,
                    probe_radius=self.vdw_probe_radius) ),
                # ( "msms_coulomb", Msms, copy_dict( msms_kwargs,
                #     probe_radius=self.probe_radius) )
            ]

            def filt( lst, lst_all ):
                if not lst:
                    return lst_all
                if lst[0]=="!":
                    return [ e for e in lst_all if e[0] not in lst[1:] ]
                else:
                    return [ e for e in lst_all if e[0] in lst ]
            self.water_variants = filt( self.variants, self.water_variants )
            self.tool_list = filt( self.tools, self.tool_list )

            for suffix, pdb_file in self.water_variants:
                for prefix, tool, tool_kwargs in self.tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ] = tool(
                        pdb_file, **copy_dict( kwargs, run=False, 
                            output_dir=self.subdir( name ), **tool_kwargs )
                    )
                    self.output_files += self.__dict__[ name ].output_files

            self.output_files += [ self.stats_file ]

    def _pre_exec( self ):
        pass
    def func( self ):
        if self.check_only:
            return

        def do( name ):
            if not self.tools:
                return True
            if self.tools[0]=="!":
                return name not in self.tools[1:]
            else:
                return name in self.tools

        if not self.analyze_only:
            if do( "opm" ):
                self.opm()
            if do( "proc" ):
                self.make_processed_pdb()
                self.make_nowat_pdb()
            if do( "dowser" ):
                self.dowser()
            if do( "msms0" ):
                self.msms0()
            if do( "dry" ):
                self.make_dry_pdb()
            if do( "final" ):
                self.make_final_pdb()
            if do( "info" ):
                self.pdb_info()
                self.opm_info()
                self.mpstruc_info()
                self.make_info()

            for suffix, pdb_file in self.water_variants:
                for prefix, tool, tool_kwargs in self.tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ]()

        # if do( "stats" ):
        #     self.make_stats()
        if do( "stats2" ):
            self.make_stats2()
        
        self.records = self.make_records()
        self.write()
        for r in self.records:
            r.info()
    def _post_exec( self ):
        if not self.parallel and not self.check_only:
            self._make_provi_file(
                pdb_id=self.pdb_id,
                pdb_file=self.relpath( self.final_pdb ),
                pdb_org_file=self.relpath( self.original_dry_pdb ),
                mplane_file=self.relpath( self.opm.mplane_file ),
                hbx_file=self.relpath( self.hbexplore_fin.hbx_file ),
                vol_file=self.relpath( self.voronoia_fin.vol_file ),
                msms_components=self.msms_vdw_fin.components_provi(
                    color="lightgreen", translucent=0.5, relpath=self.relpath ),
            )
        if self.parallel and self.check_only:
            dct = collections.defaultdict( list )
            db_records = []
            for t in self.tool_list:
                info = t.check( full=True )
                tag = re.split( "/|\.", info )[0]
                if tag!="opm":
                    # check if opm found two mplanes
                    try:
                        if len( t.opm.get_planes() )!=2:
                            tag="mplane"
                    except:
                        print tag, info, t.id
                dct[ tag ].append( t )
            for tag, t_list in dct.iteritems():
                if tag!="Ok":
                    print tag, " ".join( map( lambda x: x.id, t_list ) ) 
            for tag, t_list in dct.iteritems():
                print tag, len(t_list)
            if self.database:
                for t in dct.get( 'Ok', [] ):
                    t_info = t.get_info()
                    db_records.append( MppdDbRecord(
                        t_info['PDBID'],
                        t_info['HEADER'],
                        t_info['RESOLUTION'],
                        t_info['EXPDTA'],
                        ",".join( t_info['KEYWDS'] ),
                        t_info['OPMSpecies'],
                        t_info['OPMFamily'],
                        t_info['OPMSuperfamily'],
                        t_info['swname'],
                        t_info['swfamily'],
                        t_info['topic']
                    ))
                db = SqliteBackend( "mppd.db", MppdDbRecord )
                db.write( db_records )
            if self.extract:
                fdir = self.extract
                if not os.path.exists( fdir ):
                    os.makedirs( fdir )
                shutil.copyfile( 
                    self.outpath( "mppd.db" ),
                    os.path.join( fdir, "mppd.db" )
                )
                for t in dct.get( 'Ok', [] ):
                    flist = [
                        t.original_dry_pdb,
                        t.final_pdb,
                        t.opm.mplane_file,
                        t.hbexplore_fin.hbx_file + ".bonds",
                        t.voronoia_fin.vol_file + ".atmprop",
                        t.outpath( "mppd.provi" )
                    ]
                    flist += t.msms_vdw_fin.component_files()
                    for fsrc in flist:
                        fdst = os.path.join( fdir, t.id, t.relpath( fsrc ) )
                        if not os.path.exists( os.path.dirname( fdst ) ):
                            os.makedirs( os.path.dirname( fdst ) )
                        shutil.copyfile( fsrc, fdst )
            if self.figures:
                alpha = 0.3
                size = 7.0
                nres = collections.defaultdict( list )
                nwater = collections.defaultdict( list )
                resolution = collections.defaultdict( list )
                ncav = collections.defaultdict( list )
                sesvol = collections.defaultdict( list )
                packdens = collections.defaultdict( list )
                packdens_buried = collections.defaultdict( list )
                for t in dct.get( 'Ok', [] ):
                    stats = t.get_stats()
                    info = t.get_info()
                    for s in stats:
                        key = ( s["source"], s["segment"] )
                        nres[ key ].append( s["residue_count"] )
                        nwater[ key ].append( s["water_count"] )
                        resolution[ key ].append( 
                            try_float( info["RESOLUTION"], 0.0 )    
                        )
                        ncav[ key ].append( s["msms"] )
                        sesvol[ key ].append( s["msms_ses"] )
                        packdens[ key ].append( s["packdens_protein"] )
                        packdens_buried[ key ].append( 
                            s["packdens_protein_buried"] 
                        )
                for key in nres.keys():
                    x = np.array( nwater[ key ] )
                    y = np.array( nres[ key ] )
                    x_y = x/y
                    r = np.array( resolution[ key ] )
                    cav = np.array( ncav[ key ] )
                    cav_y = cav/y
                    vol = np.array( sesvol[ key ] ) * -1
                    vol_y = vol/y
                    pd = np.array( packdens[ key ] )
                    pd_buried = np.array( packdens_buried[ key ] )
                    
                    from mpl_toolkits.axes_grid.anchored_artists import (
                        AnchoredText
                    )

                    def hist( axis, x, label, loc=1, nzero=False ):
                        if nzero:
                            x = x[ x!=0 ]
                        if len(x)==0:
                            x = np.array([ 0 ])
                        axis.hist( x, normed=True, bins=25 )
                        axis.set_xlabel( label )
                        summary = (
                            "Var: %.4f\nStd: %.4f\nMean: %.4f\n"
                            "Median: %.4f\nMin: %.4f\nMax: %.4f\n"
                        ) % ( 
                            np.var(x), x.std(), x.mean(), 
                            np.median(x), x.min(), x.max() 
                        )
                        at = AnchoredText(
                            summary, loc=loc or 1, prop={ "size": 10 }, 
                            frameon=True, pad=0.5, borderpad=1.0
                        )
                        axis.add_artist( at )

                    def scatter( axis, x, y, xlabel, ylabel, 
                                 loc=1, nzero=True ):
                        if nzero:
                            xnzero = x!=0
                            ynzero = y!=0
                            x = x[ xnzero&ynzero ]
                            y = y[ xnzero&ynzero ]
                        try:
                            r = pearsonr(x, y)
                        except Exception:
                            r = ( np.nan, np.nan )
                        axis.scatter( x, y, alpha=alpha, s=size )
                        axis.set_xlabel( xlabel )
                        axis.set_ylabel( ylabel )
                        summary = "r: %.4f\np: %.4f\n" % r
                        at = AnchoredText(
                            summary, loc=loc or 1, prop={ "size": 10 }, 
                            frameon=True, pad=0.5, borderpad=1.0
                        )
                        axis.add_artist( at )

                    fig, (ax) = plt.subplots(3, 4, figsize=[20,12] )
                    
                    scatter( ax[0,0], x, y, "#h20", "#res" )
                    hist( ax[0,1], x_y, "#h2o / #res" )
                    
                    scatter( ax[1,0], r, x_y, "resolution [A]", "#h2o / #res" )
                    hist( ax[1,1], cav_y, "#cav / #res" )

                    hist( ax[2,0], vol_y, "ses_vol [A^3] / #res" )

                    hist( ax[0,2], pd, "packing density" )

                    scatter( ax[1,2], r, pd, 
                        "resolution [A]", "packing density" )

                    hist( ax[0,3], pd_buried, "packing density buried" )
                    
                    scatter( ax[1,3], r, pd_buried, 
                        "resolution [A]", "packing density buried" )

                    fig.savefig( "_".join( key ) + ".png" )

                def bar( ax, ydata, labels ):
                    y = [ np.array(yd).mean() for yd in ydata ]
                    x = np.arange( len( y ) )
                    e = [ np.array(yd).std() for yd in ydata ]
                    ax.bar( 
                        x, y, align='center', yerr=e, 
                        ecolor='black', facecolor='#777777' 
                    )
                    ax.set_xticks( x )
                    ax.set_xticklabels( labels )
                    xlim = ( x.min()-1, x.max()+1 )
                    ax.set_xlim( xlim )
                

                # ...
                tm_keys = [
                    ( "org", "TM" ),
                    ( "fin", "TM" ),
                    ( "dow", "TM" )
                ]
                if all( map( lambda k: k in nres, tm_keys ) ):
                    ydata = []
                    labels = []
                    for key in tm_keys:
                        ydata.append(
                            ( np.array( nwater[ key ] ) /
                                np.array( nres[ key ] ) ) * 100 
                        )
                        labels.append( key[0] )
                    fig, (ax) = plt.subplots(1, 1, figsize=[6,4.5] )
                    bar( ax, ydata, labels )
                    fig.savefig( "h2o_per_100_res.png" )


                # ...
                nwater_cutoff = collections.defaultdict( list )
                cutoff_list = np.arange(1.4, 3.9, step=0.1)
                for t in dct.get( 'Ok', [] ):
                    stats2 = t.get_stats2()
                    for s2 in stats2:
                        key = s2["segment"]
                        if not len(nwater_cutoff[ key ]):
                            for cutoff in cutoff_list:
                                nwater_cutoff[ key ].append([])
                        water_count = s2["water_count"]
                        # residue_count = s2["residue_count"]
                        count_list = s2["exp_water_cutoff_count"]
                        for i, cutoff in enumerate( cutoff_list ):
                            frac = try_div( count_list[i], water_count )
                            nwater_cutoff[ key ][i].append( frac )
                for key in nwater_cutoff.keys():
                    fig, (ax) = plt.subplots(1, 1, figsize=[8,4] )
                    bar( 
                        ax, 
                        nwater_cutoff[ key ], 
                        map( str, cutoff_list ) 
                    )
                    fig.savefig( str( key ) + ".png" )
    def make_final_pdb( self ):
        npdb_dow = NumPdb( 
            self.dowser_dry_pdb, features=self.npdb_features )
        npdb_org = NumPdb( 
            self.original_dry_pdb, features=self.npdb_features )
        dow_tree = get_tree( 
            npdb_dow.get( 'xyz', resname="HOH" ) 
        )
        sele = npdb_org.sele()
        i = 0
        for numa in npdb_org.iter_resno( incomplete=True ):
            flag = False
            if numa[0]['resname']=='HOH':
                dist = dow_tree.query( numa['xyz'][0] )[0]
                if dist>2.7:
                    flag = True
            for a in numa._atoms:
                sele[i] = flag
                i += 1
        coords_dow, atoms_dow = npdb_dow._select()
        coords_org, atoms_org = npdb_org._select( sele=sele )
        npdb_final = NumAtoms( 
            np.hstack(( atoms_dow, atoms_org )), 
            np.vstack(( coords_dow, coords_org ))
        )
        npdb_final.write2( self.final_pdb )
    def make_dry_pdb( self ):
        envelope = self.msms0.envelope_msms.get_vert( filt=False )
        envelope = [ ( x['x'], x['y'], x['z'] ) for x in envelope ]
        envelope_tree = scipy.spatial.KDTree( envelope )
        npdb_non = NumPdb( 
            self.no_water_file, features=self.npdb_features )
        non_tree = scipy.spatial.KDTree( npdb_non['xyz'] )
        
        pr2 = self.probe_radius * 2
        pr3 = self.probe_radius * 3
        wet_pdb = [
            ( self.dowser.dowser_file, self.dowser_dry_pdb ),
            ( self.processed_pdb, self.original_dry_pdb )
        ]
        for pdb_file, out_file in wet_pdb:
            npdb = NumPdb( pdb_file, features=self.npdb_features )
            sele_het = npdb.sele( record="HETATM" )
            sele_not_wat = npdb.sele( resname="HOH", invert=True )
            hetero_tree = get_tree(
                npdb.get( 'xyz', sele=( sele_het&sele_not_wat ) )
            )
            sele = npdb.sele()
            i = 0
            for numa in npdb.iter_resno( incomplete=True ):
                flag = True
                if numa[0]['resname']=='HOH':
                    # assume the first water atom is the oxygen
                    dist_env = envelope_tree.query( numa['xyz'][0] )[0]
                    dist_het = hetero_tree.query( numa['xyz'][0] )[0]
                    dist_non = non_tree.query( numa['xyz'][0] )[0]
                    # TODO het cutoff only for dowser_file
                    if dist_env<pr2 or dist_het<3.9 or dist_non>pr3:
                        flag = False
                for a in numa._atoms:
                    sele[i] = flag
                    i += 1
            npdb.copy( sele=sele ).write2( out_file )
    def make_processed_pdb( self ):
        npdb = NumPdb( 
            self.opm.processed_file, features=self.npdb_features )
        sele = npdb.sele()
        coords_dict = {}
        i = 0
        for numa in npdb.iter_resno( incomplete=True ):
            for c in numa._coords:
                c = tuple( c )
                # remove atoms with identical coords
                if c in coords_dict:
                    print npdb._atoms[i]
                    sele[i] = False
                else:
                    coords_dict[ c ] = True
                    sele[i] = True
                i += 1
        npdb.copy( sele=sele ).write2( self.processed_pdb )
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
        npdb_sol_dict = {}
        for suffix, pdb_file in self.water_variants:
            npdb = NumPdb( pdb_file, features=self.npdb_features )
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
            npdb_sol_dict[ suffix ] = npdb.copy( sele=sele )
            npdb.copy( sele=sele ).write( self.outpath( "sol_region.pdb" ) )
            np.logical_not( sele, sele )
            npdb_tm_dict[ suffix ] = npdb.copy( sele=sele )
            npdb.copy( sele=sele ).write( self.outpath( "tm_region.pdb" ) )
        return npdb_dict, npdb_tm_dict, npdb_sol_dict
    def make_stats( self ):
        segments = zip( [ "ALL", "TM", "SOL" ], self.get_npdb_dicts() )
        mplanes = self.opm.get_planes()
        dist = abs( point_plane_dist( mplanes[1][0], mplanes[0] ) )
        pr = self.probe_radius
        mppd_stats = []
        for suffix, pdb_file in self.water_variants:
            hbx = self.__dict__[ "hbexplore_%s" % suffix ]
            voro = self.__dict__[ "voronoia_%s" % suffix ]
            msms_vdw = self.__dict__[ "msms_vdw_%s" % suffix ]
            # msms_coulomb = self.__dict__[ "msms_coulomb_%s" % suffix ]
            for seg, npdb_d in segments:
                _count_water = count_water( npdb_d[ suffix ] )
                _count_hbonds = count_hbonds( npdb_d[ suffix ], hbx )
                _count_holes_voro = count_holes_voro( npdb_d[ suffix ], voro )
                _packing_density = packing_density( npdb_d[ suffix ], voro )
                _count_holes_msms = count_holes_msms( 
                    npdb_d[ suffix ], msms_vdw, pr, 
                    len( segments[0][1][ suffix ] )
                )
                mppd_stats.append({
                    "pdb_id": self.pdb_id, 
                    "source": suffix,
                    "segment": seg,
                    "mplanes_dist": dist,
                    "water_count": _count_water[0],
                    "residue_count": _count_water[1],
                    "hetero_count": _count_water[2],
                    "chain_count": _count_water[3],
                    "identical_chains": None,
                    "hbx_protein_water": _count_hbonds[0],
                    "hbx_water_water": _count_hbonds[1],
                    "hbx_water_ligand": _count_hbonds[2],
                    "voro_not_filled": _count_holes_voro[0],
                    "voro_partly_filled": _count_holes_voro[1],
                    "packdens_all": _packing_density[0],
                    "packdens_protein": _packing_density[1],
                    "packdens_water": _packing_density[2],
                    "packdens_hetero": _packing_density[3],
                    "packdens_protein_buried": _packing_density[4],
                    "msms": _count_holes_msms[0],
                    "msms_ses": _count_holes_msms[1],
                    "msms_gt_water": _count_holes_msms[2],
                    "msms_gt_water_ses": _count_holes_msms[3]
                })
        with open( self.stats_file, "w" ) as fp:
            json.dump( mppd_stats, fp )
    def make_stats2( self ):
        segments = zip( [ "ALL", "TM", "SOL" ], self.get_npdb_dicts() )
        mppd_stats2 = []
        for seg, npdb_d in segments:
            npdb_dow = npdb_d[ "dow" ]
            npdb_org = npdb_d[ "org" ]
            dow_tree = get_tree( 
                npdb_dow.get( 'xyz', resname="HOH" ) 
            )
            water_count = 0
            residue_count = 0
            cutoff_list = np.arange(1.4, 3.9, step=0.1)
            count_list = [0] * len( cutoff_list )
            for numa in npdb_org.iter_resno( incomplete=True ):
                # count org waters that are not in dow
                # at various cutoffs
                if numa[0]['resname']=='HOH':
                    water_count += 1
                    dist = dow_tree.query( numa['xyz'][0] )[0]
                    for i, cutoff in enumerate( cutoff_list ):
                        if dist>cutoff:
                            count_list[ i ] += 1
                elif numa[0]["record"]=="ATOM  ":
                    residue_count += 1
            mppd_stats2.append({
                "pdb_id": self.pdb_id, 
                "segment": seg,
                "water_count": water_count,
                "residue_count": residue_count,
                "exp_water_cutoff_count": count_list
            })
        with open( self.stats2_file, "w" ) as fp:
            json.dump( mppd_stats2, fp )
    def make_info( self ):
        info = self.pdb_dic[ self.pdb_id.upper() ]
        with open( self.info_file, "w" ) as fp:
            json.dump( info, fp, indent=4 )
    def get_info( self ):
        with open( self.info_file, "r" ) as fp:
            return json.load( fp )
    def get_stats( self ):
        with open( self.stats_file, "r" ) as fp:
            return json.load( fp )
    def get_stats2( self ):
        with open( self.stats2_file, "r" ) as fp:
            return json.load( 
                fp, object_pairs_hook=collections.OrderedDict 
            )
    def make_records( self ):
        with open( self.stats_file, "r" ) as fp:
            mppd_stats = json.load( fp )
        mppd_records = []
        for stat in mppd_stats:
            mppd_records.append( MppdRecord(
                stat["pdb_id"],
                stat["source"],
                stat["segment"],
                stat["mplanes_dist"],
                stat["water_count"],
                stat["residue_count"],
                stat["hetero_count"],
                stat["chain_count"],
                stat["identical_chains"],
                stat["hbx_protein_water"],
                stat["hbx_water_water"],
                stat["hbx_water_ligand"],
                stat["voro_not_filled"],
                stat["voro_partly_filled"],
                stat["packdens_all"],
                stat["packdens_protein"],
                stat["packdens_water"],
                stat["packdens_hetero"],
                stat["msms"],
                stat["msms_ses"],
                stat["msms_gt_water"],
                stat["msms_gt_water_ses"]
            ))
        return mppd_records


def packing_density( npdb, voronoia ):
    try:
        vol = voronoia.get_vol()
    except:
        return ( 0, 0, 0, 0, 0 )
    pd_sum = 0
    pd_sum_protein = 0
    pd_sum_protein_buried = 0
    pd_sum_water = 0
    pd_sum_hetero = 0
    count = 0
    protein_count = 0
    protein_buried_count = 0
    water_count = 0
    hetero_count = 0
    pd_dict = vol["packdens"]
    buried_dict = vol["buried"]
    for numa in npdb._iter_resno():
        for a in numa:
            pd = pd_dict.get( a["atomno"] )
            if not pd:
                continue
            pd_sum += pd
            count += 1
            if a["resname"]=="HOH":
                pd_sum_water += pd
                water_count += 1
            elif a["record"]=="ATOM  ":
                pd_sum_protein += pd
                protein_count += 1
                if buried_dict.get( a["atomno"] ):
                    pd_sum_protein_buried += pd
                    protein_buried_count += 1
            elif a["record"]=="HETATM":
                pd_sum_hetero += pd
                hetero_count += 1
    return (
        try_div( pd_sum, count ),
        try_div( pd_sum_protein, protein_count ),
        try_div( pd_sum_water, water_count ),
        try_div( pd_sum_hetero, hetero_count ),
        try_div( pd_sum_protein_buried, protein_buried_count )
    )


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
    hb_count_pw = 0
    hb_count_ww = 0
    hb_count_wl = 0
    try:
        hbonds = hbexplore.get_hbonds()
    except Exception:
        return ( 0, 0, 0 )
    for hb in hbonds:
        # check if both atoms are in the npdb
        na1 = npdb.copy( resno=hb.resno1, chain=hb.chain1 )
        na2 = npdb.copy( resno=hb.resno2, chain=hb.chain2 )
        if len(na1) and len(na2):
            if hb.type=="w-w":
                hb_count_ww += 1
            elif hb.type in ("w-l", "l-w"):
                hb_count_wl += 1
            elif hb.type in ("s-w", "w-s", "B-w", "w-B"):
                hb_count_pw += 1
    return ( hb_count_pw, hb_count_ww, hb_count_wl )


def count_holes_voro( npdb, voronoia ):
    try:
        vol = voronoia.get_vol()
    except Exception:
        return ( 0, 0 )
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


def count_holes_msms( npdb, msms, probe_radius, n ):
    try:
        area = msms.get_area()
        components = msms.get_components()
    except Exception:
        return ( 0, 0, 0, 0 )
    components0 = []
    for c in components:
        for nb_atomno in area["ses"][ c[0] ]:
            if not len( npdb.copy( atomno=nb_atomno ) ) and nb_atomno<=n:
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



class MppdStats( PyTool ):
    args = [
        _( "mppd_file", type="file", ext="json" )
    ]
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        records = JsonBackend( self.mppd_file, MppdRecord )
        print len( records )




