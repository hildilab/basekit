

from dowser import Dowser
from solvate import Solvate
from pdb import PdbInfo

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "voro_pipe" )


_VoroPipeRecord = collections.namedtuple( "_VoroPipeRecord", [
    "pdb_id", "source", "segment",
    "water_count", "residue_count", "hetero_count", 
    "chain_count",
    "voro_not_filled", "voro_partly_filled", 
    "packdens_all", "packdens_protein", "packdens_water", "packdens_hetero"
])
class VoroPipeRecord( _VoroPipeRecord ):
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



class VoroniaPipeline( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """The Voronoia pipeline"""
    args = [
        _( "pdb_input", type="str" ),
        _( "analyze_only|ao", type="checkbox", default=False ),
        _( "check_only|co", type="checkbox", default=False ),
        _( "tools", type="str", nargs="*", default=[],#tool-list= original, removed_water, dowsered, solvvated
            help="a '!' as the first arg negates the list" ),
        _( "figures|fig", type="checkbox", default=False ),
        _( "database|db", type="checkbox", default=False ),
        # voronoia shuffle
        _( "voro_shuffle", type="checkbox", default=False ),
        # dowser max repeats
        _( "dowser_max", type="int", default=None ),
    ]
    out = [
        _( "processed_pdb", file="proc.pdb" ),
        _( "no_water_file", file="nowat.pdb" ),
        _( "dowser_pdb", file="dowser.pdb" ),
        _( "final_pdb", file="final.pdb" ),
        _( "stats_file", file="stats.json" ),
        _( "stats2_file", file="stats2.json" ),
        _( "info_file", file="info.json" ),
    ]
    RecordsClass = VoroPipeRecord
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

            self.pdb_info = PdbInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("pdb_info") ) )
            self.output_files += [ self.pdb_info.info_file ]
            
            self.splitted = Splitted(
                    self.pdb_id,
                    **copy_dict( kwargs, run=False, 
                        output_dir=self.subdir("splitted") )
                )
            
            self.moderna = Moderna(
                    self.pdb_id,
                    **copy_dict( kwargs, run=False, 
                        output_dir=self.subdir("moderna") )
                )
            
            if self.use_rem_wat:
                self.dry = Remove_Water(############################
                    self.pdb_id,
                    **copy_dict( kwargs, run=False, 
                        output_dir=self.subdir("removed_water") )
                )
            elif self.use_dowser:
                self.dowser = DowserRepeat(
                self.pdb_id,
                **copy_dict( kwargs, run=False, alt='x',
                    output_dir=self.subdir("dowser"),
                    max_repeats=self.dowser_max )
            )
            elif self.use_solvate:
                self.solvate = Solvate(
                    self.pdb_id,
                    **copy_dict( kwargs, run=False, 
                        output_dir=self.subdir("solvated") )
                )
            else:
                self.original = Copy_to_dir(###########################
                    self.pdb_id,
                    **copy_dict( kwargs, run=False, 
                        output_dir=self.subdir("original") )
                )
            
            self.output_files += [ self.dry.output_files, self.dowser.output_files ]
            self.output_files += [ self.solvate.output_files, self.original.output_files ]

            msms_kwargs = { 
                "all_components": True,
                "density": 1.0, "hdensity": 3.0,
                "envelope": self.probe_radius * 2,
                "envelope_hclust": self.envelope_hclust,
                "atom_radius_add": self.atom_radius_add,
            }
            self.msms0 = Msms(
                self.no_water_file, **copy_dict( kwargs, run=False, 
                output_dir=self.subdir( "msms0" ), 
                **copy_dict( msms_kwargs, probe_radius=self.vdw_probe_radius,
                    all_components=False )
            ))
            self.output_files += self.msms0.output_files
            self.output_files += [ self.original_dry_pdb, self.final_pdb ]

            self.dssp = Dssp(
                self.final_pdb, **copy_dict( kwargs, run=False, 
                output_dir=self.subdir( "dssp" )
            ))
            self.output_files += self.dssp.output_files

            self.opm_info = OpmInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("opm_info") ) )
            self.mpstruc_info = MpstrucInfo( self.pdb_id,
                **copy_dict( kwargs, run=False, 
                    output_dir=self.subdir("mpstruc_info") ) )

            self.water_variants = [
                ( "non", self.no_water_file ),
                ( "org", self.original_dry_pdb ),
                ( "dow", self.dowser_dry_pdb ),
                ( "fin", self.final_pdb ),
            ]
            
            self.tool_list = [
                ( "voronoia", Voronoia, { 
                    "ex":0.2, "shuffle": self.voro_shuffle
                }),
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
            self.do_tool_list = filt( self.tools, self.tool_list )

            for suffix, pdb_file in self.water_variants:
                for prefix, tool, tool_kwargs in self.tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ] = tool(
                        pdb_file, **copy_dict( kwargs, run=False, 
                            output_dir=self.subdir( name ), **tool_kwargs )
                    )
                    self.output_files += self.__dict__[ name ].output_files

            
            self.output_files += [ self.mpstruc_info.info_file ]
            self.output_files += [ self.opm_info.info_file ]
            self.output_files += [ self.info_file ]
            self.output_files += [ self.outpath( "mppd.provi" ) ]
            self.output_files += [ self.stats_file ]
            self.output_files += [ self.stats2_file ]

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
###prepair PDB for calculation###
#get_pdb
#splitted, moderna
#original PDB
#removed water pdb
#dowsered pdb
#solvated pdb
###calculate get_volume###
###do Analysis###
        if not self.analyze_only:
            if do( "get_pdb" ):
                self.get_pdb()
            if do( "splitted" ):
                self.splitted()
            if do( "moderna" ):
                self.moderna()
            if do( "dry" ):
                self.make_dry_pdb()
            if do( "dowsered_pdb" ):
                self.dowsered_pdb()
            if do( "solvated_pdb" ):
                self.solvated_pdb()
            if do( "pdb_info" ):
                self.pdb_info()
            if do( "info" ):
                self.make_info()
            if do( "provi" ):
                npdb = NumPdb( self.final_pdb, features=self.npdb_features )
                self._make_provi_file(
                    pdb_id=self.pdb_id,
                    pdb_title=self.pdb_info.get_info()["title"],
                    pdb_file=self.relpath( self.final_pdb ),
                    pdb_org_file=self.relpath( self.original_dry_pdb ),
                    mplane_file=self.relpath( self.opm.mplane_file ),
                    hbx_file=self.relpath( self.hbexplore_fin.hbx_file ),
                    vol_file=self.relpath( self.voronoia_fin.vol_file ),
                    msms_components=self.msms_vdw_fin.components_provi(
                        color="lightgreen", translucent=0.5, 
                        relpath=self.relpath, max_atomno=len( npdb )
                    ),
                )

            for suffix, pdb_file in self.water_variants:
                for prefix, tool, tool_kwargs in self.do_tool_list:
                    name = "%s_%s" % ( prefix, suffix )
                    self.__dict__[ name ]()

        if do( "stats" ):
            self.make_stats()
        if do( "stats2" ):
            self.make_stats2()
        
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
                if tag=="mppd":
                    tag = "provi"
                if tag!="opm":
                    # check if opm found two mplanes
                    try:
                        if len( t.opm.get_planes() )!=2:
                            tag = "mplane"
                    except:
                        #print tag, check_info, t.id
                        pass
                else:
                    if os.path.isfile( t.opm.outpath( "ppm_error.txt" ) ):
                        tag = "ppm"
                        # print open( t.opm.outpath( "ppm_error.txt" ) ).read()
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


            # representative id search
            test_rep_list = flatten([
                zip( [x]*len(dct[x]), dct[x] ) 
                for x in [ "opm", "ppm", "msms0", "msms_vdw_fin", "dowser" ]
            ])
            for tag, t in test_rep_list:
                opm_info = t.opm_info.get_info()
                mpstruc_info = t.mpstruc_info.get_info()
                rid_list = []
                if opm_info:
                    rep_id = opm_info.get("representative")
                    if rep_id:
                        rid_list.append( rep_id.upper() )
                    rid_list += opm_info.get("related_ids", [])
                if mpstruc_info:
                    master_id = mpstruc_info.get("master")
                    if master_id:
                        rid_list.append( master_id.upper() )
                    rid_list += mpstruc_info.get("related", [])
                        
                rep = None
                for rid in rid_list:
                    for x in dct["Ok"]:
                        if x.pdb_id==rid:
                            rep = x
                            break
                    else:
                        continue
                    break
                if rep:
                    dct[ "representative" ].append( t )
                else:
                    dct[ "no_representative" ].append( t )

            ignore_tags = [ 
                "no_pdb_entry",
                "mplane",
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
                elif tag in [ "opm", "ppm", "msms0", "msms_vdw_fin" ]:
                    if pdb_id in dct[ "representative" ]:
                        status = "linked"
                    else:
                        status = "pending"
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

                            t_info['opm_superfamily'],
                            t_info['opm_family'],
                            t_info['opm_representative'],
                            t_info['opm_species'],
                            ",".join( t_info['opm_related'] ),

                            t_info['mpstruc_group'],
                            t_info['mpstruc_subgroup'],
                            t_info['mpstruc_name'],
                            t_info['mpstruc_species'],
                            t_info['mpstruc_master'],
                            ",".join( t_info['mpstruc_related'] ),

                            cur_info.get('representative', ""),
                            ",".join( cur_info.get('related', []) ),

                            status_dct[ t.pdb_id ],

                            t_stats.get('packdens_protein_buried'),
                            t_stats.get('water_count'),
                            t_stats.get('residue_count'),
                            t_stats.get('msms'),
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
                        if s["segment"]!="TM":
                             continue
                        key = ( s["source"], s["segment"] )
                        nres[ key ].append( s["residue_count"] )
                        nwater[ key ].append( s["water_count"] )
                        resolution[ key ].append( 
                            try_float( info["pdb_resolution"], 0.0 )    
                        )
                        ncav[ key ].append( s["msms"] )
                        sesvol[ key ].append( s["msms_ses"] )
                        packdens[ key ].append( s["packdens_protein"] )
                        packdens_buried[ key ].append( 
                            s["packdens_protein_buried"] 
                        )
                print nres.keys()
                for key in nres.keys():
                    print key
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
                        axis.set_ylim( ( 0, axis.get_ylim()[1] ) )
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
    def make_nowat_pdb( self ):
        with open( self.no_water_file, "w" ) as fp:
            with open( self.processed_pdb, "r" ) as fp_pdb:
                for line in fp_pdb:
                    if ( line[0:6] in ["ATOM  ", "HETATM"] and
                            line[17:20]!="HOH" ):
                        fp.write( line )
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
        pdb_info = self.pdb_info.get_info() or {}
        opm_info = self.opm_info.get_info() or {}
        mpstruc_info = self.mpstruc_info.get_info() or {}
        info = {
            "pdb_id": self.pdb_id.upper(),
            "pdb_title": pdb_info.get("title"),
            "pdb_keywords": pdb_info.get("keywords", []),
            "pdb_experiment": pdb_info.get("experiment"),
            "pdb_resolution": pdb_info.get("resolution"),

            "opm_superfamily": opm_info.get("superfamily"),
            "opm_family": opm_info.get("family"),
            "opm_representative": opm_info.get("representative"),
            "opm_species": opm_info.get("species"),
            "opm_related": opm_info.get("related_ids", []),

            "mpstruc_group": mpstruc_info.get("group"),
            "mpstruc_subgroup": mpstruc_info.get("subgroup"),
            "mpstruc_name": mpstruc_info.get("name"),
            "mpstruc_species": mpstruc_info.get("species"),
            "mpstruc_master": mpstruc_info.get("master"),
            "mpstruc_related": mpstruc_info.get("related", []),
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


        