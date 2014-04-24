from __future__ import with_statement
from __future__ import division


import os
import re
import math
import json
import shutil

import collections
from dowser import DowserRepeat
from solvate import Solvate
from moderna_tools import Moderna
from pdb import PdbInfo
from pdb import JoinSplitted
from voronoia import Voronoia
from pdb import get_pdb_files
from utils.tool import (
    _, _dir_init, PyTool, RecordsMixin, ParallelMixin, ProviMixin,
    JsonBackend, SqliteBackend
)
from utils import copy_dict, try_float, try_div, flatten, dir_walker
from utils.numpdb import NumPdb, NumAtoms

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "voro_pipe" )


_VoroRecord = collections.namedtuple( "_VoroPipeRecord", [
    "pdb_id", "source", "segment",
    "water_count", "residue_count", "hetero_count", 
    "chain_count",
    "voro_not_filled", "voro_partly_filled", 
    "packdens_all", "packdens_protein", "packdens_water", "packdens_hetero"
])
class VoroRecord( _VoroRecord ):
    def info( self ):
        print "### %s [ %s, %s ]" % ( 
            self.pdb_id, self.source, self.segment 
        )
        print "counts: %i (water), %i (aa), %i (hetero)" % (
            self.water_count, self.residue_count, self.hetero_count
        )
        print "chains: %i" % (
            self.chain_count
        )
        print "voro: %i (not filled), %i (partly filled)" % (
            self.voro_not_filled, self.voro_partly_filled
        )
        print "pd: %0.3f (all), %0.3f (aa), %0.3f (water), %0.3f (hetero)" % (
            self.packdens_all, self.packdens_protein, 
            self.packdens_water, self.packdens_hetero
        )
        print ""



class VoronoiaPipeline( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """The Voronoia pipeline"""
    args = [
        _( "pdb_input", type="str" ),
        _( "analyze_only|ao", type="bool", default=False ),
        _( "check_only|co", type="bool", default=False ),
        _( "variants", type="str", nargs="*", default=[], 
            help="a '!' as the first arg negates the list" ),
        _( "tools", type="str", nargs="*", default=[],
            help="a '!' as the first arg negates the list" ),
        _( "figures|fig", type="bool", default=False ),
        _( "database|db", type="bool", default=False ),
        # voronoia shuffle
        _( "voro_shuffle", type="bool", default=False ),
        # dowser max repeats
        _( "dowser_max", type="int", default=None ),
    ]
    out = [
        _( "moderna_file", file="moderna.pdb", optional=True  ),
       # _( "no_water_file", file="nowat.pdb", optional=True  ),
        _( "dowser_file", file="dowser.pdb", optional=True  ),
        _( "original_file", file="original.pdb", optional=True  ),
       # _( "final_file", file="final.pdb", optional=True  ),
        _( "solvate_file", file="solvate.pdb", optional=True  ),
        _( "pd_at_res", file="pd_at_res.json", optional=True  ),
       # _( "stats_file", file="stats.json", optional=True  ),
       # _( "stats2_file", file="stats2.json", optional=True  ),
       # _( "info_file", file="info.json", optional=True  ),
    ]
    RecordsClass = VoroRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "voro.provi"
    def _init( self, *args, **kwargs ):
        self._init_records( None, **kwargs )
        self._init_parallel( self.pdb_input, **kwargs )

        if not self.parallel:
            self.pdb_id = self.pdb_input
            self.npdb_features = {
                "sstruc": False,
                "phi_psi": False
            }

            self.output_files = []
            shutil.copyfile(self.pdb_input, self.original_file)
            self.pdb_info = PdbInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("pdb_info") ) )
            self.output_files += [ self.pdb_info.info_file ]
            
            #if self.pdb_info.get_info()["splited_entry"]:
            #    self.splitted = JoinSplitted(
            #            self.pdb_info["splited_entry"],
            #            **copy_dict( kwargs, run=False, 
            #                output_dir=self.subdir("splitted") )
            #        )
            #    self.output_files += self.splitted.output_files
            
            self.moderna = Moderna(
                    self.pdb_id, ws=True, tool="clean_structures",
                    **copy_dict( kwargs, run=False, 
                        output_dir=self.subdir("moderna") )
                )
            self.output_files += self.moderna.output_files
            for elem in self.moderna.output_files:
                if "moderna.pdb" in elem:
                    self.moderna_file = elem
            
            self.dowser = DowserRepeat(
                self.moderna_file,
                **copy_dict( kwargs, run=False, alt='x',
                    output_dir=self.subdir("dowser"),
                    max_repeats=self.dowser_max )
            )
            self.output_files += self.dowser.output_files
            for elem in self.dowser.output_files:
                if "dowser.pdb" in elem:
                    self.dowser_file = elem
            
            self.solvate = Solvate(
                self.moderna_file,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("solvate") )
            )
            self.output_files += self.solvate.output_files
            for elem in self.solvate.output_files:
                if "solvate.pdb" in elem:
                    self.solvate_file = elem
            self.output_files += self.original_file
            
            self.pdb_variants = [
                ( "mod", self.moderna_file ),
                ( "org", self.original_file ),
                ( "dow", self.dowser_file ),
                ( "sol", self.solvate_file )
            ]
            
            
            self.tool_list = [
                ( "voronoia", Voronoia, { 
                    "ex":0.2, "shuffle": self.voro_shuffle,
                    "get_nrholes":True, "make_reference":True
                })
            ]

            def filt( lst, lst_all ):
                if not lst:
                    return lst_all
                if lst[0]=="!":
                    return [ e for e in lst_all if e[0] not in lst[1:] ]
                else:
                    return [ e for e in lst_all if e[0] in lst ]
                
            self.pdb_variants = filt( self.variants, self.pdb_variants )
            self.do_tool_list = filt( self.tools, self.tool_list )

            for suffix, pdb_file in self.pdb_variants:
                for prefix, tool, tool_kwargs in self.tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ] = tool(
                        pdb_file, **copy_dict( kwargs, run=False, 
                            output_dir=self.subdir( name ), **tool_kwargs )
                    )
                    self.output_files += self.__dict__[ name ].output_files

            
            #self.output_files += [ self.info_file ]
            #self.output_files += [ self.outpath( "voro.provi" ) ]
            #self.output_files += [ self.stats_file ]
            #self.output_files += [ self.stats2_file ]

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
            #if do( "get_pdb" ):
            #    self.get_pdb_file()
            #if do( "splitted" ):
            #    self.splitted()
            if do( "moderna" ):
                self.moderna()
            if do( "dowser" ):
                self.dowser()
            if do( "solvate" ):
                self.solvate()
            if do( "final" ):
                self.clean_up_pdbs()
        #    if do( "final" ):
        #        self.make_final_pdb()
        #    if do( "pdb_info" ):
        #        self.pdb_info()
        #    if do( "info" ):
        #        self.make_info()
        #    if do( "provi" ):
        #        npdb = NumPdb( self.final_file, features=self.npdb_features )
        #        self._make_provi_file(
        #            pdb_id=self.pdb_id,
        #            pdb_title=self.pdb_info.get_info()["title"],
        #            pdb_file=self.relpath( self.final_file ),
        #            pdb_org_file=self.relpath( self.original_file ),
        #            vol_file=self.relpath( self.voronoia_fin.vol_file ),
        #        )
            for suffix, pdb_file in self.pdb_variants:
                for prefix, tool, tool_kwargs in self.do_tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ]()
        #
        #if do( "stats" ):
        #    self.make_stats()
        #if do( "stats2" ):
        #    self.make_stats2()
        
        #self.records = self.make_records()
        self.write()
        # for r in self.records:
        #     r.info()
    def _post_exec( self ):
        def merge(obj_1, obj_2):
            if type(obj_1) == dict and type(obj_2) == dict:
                result = {}
                for key, value in obj_1.iteritems():
                    if key not in obj_2:
                        result[key] = value
                    else:
                        result[key] = merge(value, obj_2[key])
                for key, value in obj_2.iteritems():
                    if key not in obj_1:
                        result[key] = value
                return result
            if type(obj_1) == list and type(obj_2) == list:
                return obj_1 + obj_2
            return obj_2
        if self.parallel and self.check_only:
            dct = collections.defaultdict( list )
            status_dct = {}
            tag_dct = {}
            db_records = []
            for t in self.tool_list:
                check_info = t.check( full=True )
                cur_info = CURATED_INFO.get( t.pdb_id, {} )
                tag = re.split( "/|\.", check_info )[0]
                if cur_info.get("no_pdb_entry"):
                    tag = "no_pdb_entry"
                tag_dct[ t.pdb_id ] = tag
                dct[ tag ].append( t )


            # representative id search
            test_rep_list = flatten([
                zip( [x]*len(dct[x]), dct[x] ) 
                for x in [ "dowser" ]
            ])
            ignore_tags = [ 
                "no_pdb_entry"
            ]

            # status types
            #   included:   all good
            #   linked:     only a representative available
            #   pending:    to be included
            #   obsolete:   superseeded
            #   defect:     low resolution; missing atoms
            #   model:      theoretical model
            for pdb_id, tag in tag_dct.iteritems():
                if tag=="Ok":
                    status = "included"
                elif tag in ignore_tags:
                    continue
                else:
                    status = "unknown"
                    print tag
                status_dct[ pdb_id ] = status


            for tag, t_list in dct.iteritems():
                if tag!="Ok":
                    print tag, " ".join( map( lambda x: x.id, t_list ) ) 
            for tag, t_list in dct.iteritems():
                print tag, len(t_list)

            if self.database:
                for tag, t_list in dct.iteritems():
                    if tag in ignore_tags:
                        continue
                    for t in t_list:
                        t_info = t.get_info()
                        if tag=="Ok":
                            for s in t.get_stats():
                                if s['source']=='fin' and s['segment']=='TM':
                                    t_stats = s
                                    break
                            else:
                                raise Exception('no stats found %s' % t.pdb_id)
                        else:
                            t_stats = {}
                        cur_info = CURATED_INFO.get( t.pdb_id, {} )
                        db_records.append( MppdDbRecord(
                            t_info['pdb_id'],
                            t_info['pdb_title'],
                            ",".join( t_info['pdb_keywords'] ),
                            t_info['pdb_experiment'],
                            t_info['pdb_resolution'],

                            status_dct[ t.pdb_id ],

                            t_stats.get('packdens_protein_buried'),
                            t_stats.get('water_count'),
                            t_stats.get('residue_count'),
                        ))
                db = SqliteBackend( "voro.db", VoroDbRecord )
                db.write( db_records )
            if self.extract:
                fdir = self.extract
                if not os.path.exists( fdir ):
                    os.makedirs( fdir )
                shutil.copyfile( 
                    self.outpath( "voro.db" ),
                    os.path.join( fdir, "voro.db" )
                )
                for t in dct.get( 'Ok', [] ):
                    flist = [
                        t.original_pdb,
                        t.final_pdb,
                        t.voronoia_fin.vol_file + ".atmprop",
                        t.outpath( "voro.provi" )
                    ]
                    flist += t.msms_vdw_fin.component_files()
                    for fsrc in flist:
                        fdst = os.path.join( fdir, t.id, t.relpath( fsrc ) )
                        if not os.path.exists( os.path.dirname( fdst ) ):
                            os.makedirs( os.path.dirname( fdst ) )
                        shutil.copyfile( fsrc, fdst )
        if self.parallel:
            p = "pd_at_res.json"
            dic_ref={}
            for m, filepath in dir_walker( self.output_dir, p ):
                with open (filepath, 'r') as fp:
                    fi= json.load(fp)
                    dic_ref = merge(dic_ref, fi)
            for key1 in dic_ref:
                for key2 in dic_ref[key1]:
                    dic_ref[key1][key2]=sum(dic_ref[key1][key2])/len(dic_ref[key1][key2])
            with open( self.pd_at_res, "w" ) as fp:
                json.dump( dic_ref, fp, indent=4 )
    def clean_up_pdbs(self):
        for old_file in [self.moderna_file, self.dowser_file, self.solvate_file]:
            temp_file=self.original_file+'_tmp.pdb'
            with open(temp_file, 'w') as fp2:
                with open(self.original_file, 'r') as fp:
                    for line in fp:
                        if line.startswith('ATOM'):
                            break
                        else:
                            fp2.write(line)
                with open(old_file, 'r') as fp:
                    for line in fp:
                        if line.startswith('REMARK'):
                            pass
                        else:
                            fp2.write(line)
            shutil.copyfile( temp_file, old_file )
            os.remove(temp_file)
    
    def make_final_pdb( self ):
        npdb_dow = NumPdb( 
            self.dowser_file, features=self.npdb_features )
        npdb_org = NumPdb( 
            self.original_file, features=self.npdb_features )
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
        npdb_final.write2( self.final_file )

    def make_processed_pdb( self ):
        npdb = NumPdb( 
            self.moderna.moderna_file, features=self.npdb_features )
        sele = npdb.sele()
        coords_dict = {}
        i = 0
        tree = scipy.spatial.KDTree( npdb['xyz'] )
        for numa in npdb._iter_resno():
            for c in numa._coords:
                c = tuple( c )
                # remove atoms with identical coords
                if c in coords_dict:
                    print npdb._atoms[i]
                    sele[i] = False
                else:
                    coords_dict[ c ] = True
                    sele[i] = True
                    # remove atoms with almost identical coords
                    rslt = tree.query( c, k=10, distance_upper_bound=0.2 )
                    rslt[1].sort()
                    if rslt[1][0]!=i and rslt[0][1]!=np.inf:
                        print i, npdb._atoms[i]
                        sele[i] = False
                i += 1
        npdb.copy( sele=sele ).write2( self.processed_pdb )

    def get_npdb_dicts( self ):
        mplanes = np.array( self.opm.get_planes() )
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
        mplanes = np.array( self.opm.get_planes() )
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

                    "water_count": _count_water[0],
                    "residue_count": _count_water[1],
                    "hetero_count": _count_water[2],
                    "chain_count": _count_water[3],
                    "identical_chains": None,

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
        pdb_info = self.pdb_info.get_info() or {}
        opm_info = self.opm_info.get_info() or {}
        mpstruc_info = self.mpstruc_info.get_info() or {}
        info = {
            "pdb_id": self.pdb_id.upper(),
            "pdb_title": pdb_info.get("title"),
            "pdb_keywords": pdb_info.get("keywords", []),
            "pdb_experiment": pdb_info.get("experiment"),
            "pdb_resolution": pdb_info.get("resolution"),
        }
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
            voro_stats = json.load( fp )
        voro_records = []
        for stat in voro_stats:
            voro_records.append( VoroRecord(
                stat["pdb_id"],
                stat["source"],
                stat["segment"],
                stat["water_count"],
                stat["residue_count"],
                stat["hetero_count"],
                stat["chain_count"],
                stat["identical_chains"],
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
        return voro_records


        