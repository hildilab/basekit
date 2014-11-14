from __future__ import with_statement
from __future__ import division

import os
import re
import json
import collections

import numpy as np
import scipy.spatial

from utils import copy_dict, try_div
from utils.tool import (
    _, _dir_init, PyTool, RecordsMixin, ParallelMixin, ProviMixin,
    JsonBackend, SqliteBackend
)
from utils.numpdb import NumPdb
from voronoia import Voronoia, HOLE_NOT_FILLED, HOLE_PARTLY_FILLED
from pdb import PdbInfo


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "voronoia" )


CURATED_INFO = collections.defaultdict( dict, {
    "": {
        "backbone_only": True,
        "no_pdb_entry": True,
        "representative": "",
        "related": [],
        "species": "",
        "comment": ""
    },
    "2XND": {
        "representative": "",
        # from sequence similarity
        "related": [ "2XOK", "2WPD", "3ZRY", "4B2Q" ],
        "comment": "membrane plane by superposition with 2WPD"
    },
    "4FZ0": {
        "comment": "membrane plane by superposition with 2QTS"
    },
    "ROGL": {
        "no_pdb_entry": True,
        "comment": ""
    },
})


for x in [
    # outer membrane complex of a type IV secretion system
    "3JQO",
]:
    CURATED_INFO[ x ].update({
        "no_transmembrane": True
    })


_VoronoiaRecord = collections.namedtuple( "_VoronoiaRecord", [
    "pdb_id", "source", "segment", 
    "water_count", "residue_count", "hetero_count", 
    "chain_count", "identical_chains",
    "voro_not_filled", "voro_partly_filled", 
    "packdens_all", "packdens_protein", "packdens_water", "packdens_hetero",
    "pdb_zscorerms"
])
class VoronoiaRecord( _VoronoiaRecord ):
    def info( self ):
        print "### %s [ %s, %s ]" % ( 
            self.pdb_id, self.source, self.segment 
        )
        print "counts: %i (water), %i (aa), %i (hetero)" % (
            self.water_count, self.residue_count, self.hetero_count
        )
        print "chains: %i, identical: %s" % (
            self.chain_count, "??? (TODO)"
        )
        print "voro: %i (not filled), %i (partly filled)" % (
            self.voro_not_filled, self.voro_partly_filled
        )
        print "pd: %0.3f (all), %0.3f (aa), %0.3f (water), %0.3f (hetero)" % (
            self.packdens_all, self.packdens_protein, 
            self.packdens_water, self.packdens_hetero
        )
        print "voro: %i" % (
            self.pdb_zscorerms
        )
        print ""


VoronoiaDbRecord = collections.namedtuple( "VoronoiaDbRecord", [
    "pdb_id", "pdb_title", "pdb_keywords", "pdb_experiment",
   # "pdb_header",
    "pdb_resolution",
    "pdb_zscorerms",
    "curated_representative", "curated_related",
    "status",
    "tm_packdens_protein_buried", "tm_water_count", "tm_residue_count", 
])



def get_tree( coords ):
    if len( coords )==0:
        coords = np.array([[ np.inf, np.inf, np.inf ]])
    return scipy.spatial.KDTree( coords )


class VoronoiaPipeline( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """The Voronoia pipeline"""
    args = [
        _( "pdb_input", type="str" ),
        _( "ex", type="float", range=[0.01, 0.5], step=0.01, default=0.1 ),
        _( "analyze_only|ao", type="bool", default=False ),
        _( "check_only|co", type="bool", default=False ),
        _( "tools", type="str", nargs="*", default=[],
            help="a '!' as the first arg negates the list" ),
        _( "database|db", type="bool", default=False ),
        # voronoia shuffle
        _( "voro_shuffle", type="bool", default=False ),
    ]
    out = [
        _( "stats_file", file="stats.json", optional=True ),
        _( "info_file", file="info.json", optional=True ),
    ]
    RecordsClass = VoronoiaRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "voronoia.provi"
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
            self.pdb_info = PdbInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("pdb_info") ) )
            self.output_files += [ self.pdb_info.info_file ]

            self.water_variants = [
                ( "ori", self.pdb_input),
            ]
            
            self.tool_list = [
                ( "voronoia", Voronoia, { 
                    "ex":self.ex, "shuffle": self.voro_shuffle,
                    "make_reference":True, "get_nrholes":True
                }),
            ]

            def filt( lst, lst_all ):
                if not lst:
                    return lst_all
                if lst[0]=="!":
                    return [ e for e in lst_all if e[0] not in lst[1:] ]
                else:
                    return [ e for e in lst_all if e[0] in lst ]
            
            self.do_tool_list = filt( self.tools, self.tool_list )
            for suffix, pdb_file in self.water_variants:
                for prefix, tool, tool_kwargs in self.tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ] = tool(
                        pdb_file, **copy_dict( kwargs, run=False, 
                            output_dir=self.subdir( name ), **tool_kwargs )
                    )
                    self.output_files += self.__dict__[ name ].output_files

            self.output_files += [ self.info_file ]
            self.output_files += [ self.outpath( "voronoia.provi" ) ]
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
            if do( "pdb_info" ):
                self.pdb_info()
            if do( "info" ):
                self.make_info()
            if do( "provi" ):
                npdb = NumPdb( self.pdb_input, features=self.npdb_features )
                self._make_provi_file(
                    pdb_id=self.pdb_id,
                    pdb_title=self.pdb_info.get_info()["title"],
                    pdb_file=self.relpath( self.pdb_input ),
                    vol_file=self.relpath( self.voronoia_ori.vol_file ),
                )

            for suffix, pdb_file in self.water_variants:
                for prefix, tool, tool_kwargs in self.do_tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ]()

        if do( "stats" ):
            self.make_stats()
            
        self.records = self.make_records()
        self.write()
        # for r in self.records:
        #     r.info()
    def _post_exec( self ):
        if self.parallel and self.check_only:
            dct = collections.defaultdict( list )
            status_dct = {}
            tag_dct = {}
            db_records = []
            for t in self.tool_list:
                check_info = t.check( full=True )
                cur_info = CURATED_INFO.get( t.pdb_id, {} )
                tag = re.split( "/|\.", check_info )[0]

                if tag=="voronoia":
                    tag = "provi"

                if tag!="pdb_info" and t.pdb_info.check():
                    info = t.pdb_info.get_info()
                    if "CA ATOMS ONLY" in info.get( "model_type", {} ):
                        tag = "calpha_only"
                    if cur_info.get("backbone_only"):
                        tag = "backbone_only"
                    if "THEORETICAL MODEL"==info.get( "experiment", "" ):
                        tag = "theoretical_model"
                    res = info.get( "resolution" )
                    if res and res>=4.0 and res!="NOT":
                        tag = "resolution"
                    if info.get( "obsolete" ):
                        tag = "obsolete"

                if cur_info.get("no_pdb_entry"):
                    tag = "no_pdb_entry"
                    
                if cur_info.get("no_transmembrane"):
                    tag = "no_transmembrane"

                tag_dct[ t.pdb_id ] = tag
                dct[ tag ].append( t )

            ignore_tags = [ 
                "no_pdb_entry",
                "representative",
                "no_representative",
                "no_transmembrane",
                "theoretical_model",
                "obsolete",
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
                elif tag=="obsolete":
                    status = "obsolete"
                elif tag=="theoretical_model":
                    status = "model"
                elif tag in [ "calpha_only", "backbone_only" ]:
                    status = "backbone_only"
                elif tag=="resolution":
                    status = "low_resolution"
                elif tag in ignore_tags:
                    continue
                else:
                    status = "unknown"
                    print tag
                status_dct[ pdb_id ] = status

            for tag, t_list in dct.iteritems():
                if tag!="Ok":
                    #include new entrys
                    #missing = map( lambda x: x.id, t_list )
                    #for elem in missing:
                    #    print elem
                    #    voropipe = VoronoiaPipeline(
                    #        os.path.join(self.pdb_input, elem + '.pdb'),
                    #        run=False, 
                    #        output_dir=os.path.join(self.output_dir, 'parallel', elem)
                    #    )
                    #    voropipe()
                    missings = map( lambda x: x.id, t_list )
                    missing = map( lambda x: os.path.join(self.pdb_input, x + '.pdb'), missings)
                    #print missing
                    for elem in missing:
                        print elem
                        voropipe = VoronoiaPipeline(
                            elem,
                    #        str(missing),
                    #        parallel='list',
                            database=True,
                            ex=self.ex,
                            run=False,
                            verbose=False,#self.verbose,
                            output_dir=os.path.join(self.output_dir, 'parallel', elem[-8:-4])
                            #output_dir=self.output_dir#os.path.join(self.output_dir, 'parallel', elem)
                        )
                        voropipe()
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
                                if s['source']=='ori' and s['segment']=='TM':
                                    t_stats = s
                                    break
                            else:
                                raise Exception('no stats found %s' % t.pdb_id)
                        else:
                            t_stats = {}
                        cur_info = CURATED_INFO.get( t.pdb_id, {} )
                        db_records.append( VoronoiaDbRecord(
                            t_info['pdb_id'],
                            t_info['pdb_title'],
                            ",".join( t_info['pdb_keywords'] ),
                            t_info['pdb_experiment'],
                            t_info['pdb_resolution'],
                            t_stats.get('pdb_zscorerms'),
                            cur_info.get('representative', ""),
                            ",".join( cur_info.get('related', []) ),
                            status_dct[ t.pdb_id ],
                            t_stats.get('packdens_protein_buried'),
                            t_stats.get('water_count'),
                            t_stats.get('residue_count'),
                        ))
                db = SqliteBackend( "voronoia.sqlite", VoronoiaDbRecord )
                db.write( db_records )

    def get_npdb_dicts( self ):
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
        voronoia_stats = []
        for suffix, pdb_file in self.water_variants:
            voro = self.__dict__[ "voronoia_%s" % suffix ]
            for seg, npdb_d in segments:
                _pdb_zscorerms = 0#pdb_zscorerms( npdb_d[ suffix ], voro )
                _count_water = count_water( npdb_d[ suffix ] )
                _count_holes_voro = count_holes_voro( npdb_d[ suffix ], voro )
                _packing_density = packing_density( npdb_d[ suffix ], voro )
                _pdb_zscorerms = self.get_voro_info()
                voronoia_stats.append({
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
                    "pdb_zscorerms": _pdb_zscorerms["pdb_zscorerms"],
                })
        with open( self.stats_file, "w" ) as fp:
            json.dump( voronoia_stats, fp )

    def make_info( self ):
        pdb_info = self.pdb_info.get_info() or {}
        info = {
            "pdb_id": self.pdb_id.upper().split("/")[-1][:-4],
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

    def get_voro_info( self ):
        for index, elem in enumerate(self.output_files):
            if "voronoia_records.json" in elem:
                with open( self.output_files[index], "r" ) as fp:
                    return json.load( fp )[0]

    def make_records( self ):
        with open( self.stats_file, "r" ) as fp:
            voronoia_stats = json.load( fp )
        voronoia_records = []
        for stat in voronoia_stats:
            voronoia_records.append( VoronoiaRecord(
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
                stat["pdb_zscorerms"]
            ))
        return voronoia_records


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


class VoronoiaStats( PyTool ):
    args = [
        _( "voronoia_file", type="file", ext="json" )
    ]
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        records = JsonBackend( self.voronoia_file, VoronoiaRecord )
        print len( records )




