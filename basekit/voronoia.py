from __future__ import with_statement
from __future__ import division



import os
import re
import json
import logging
import collections
import utils.numpdb as numpdb
import math

import numpy as np
from numpy import *
from array import array

from basekit import utils
from utils import memoize_m
from utils.tool import _, _dir_init, CmdTool, ProviMixin, ParallelMixin
from utils.tool import RecordsMixin, PyTool

import provi_prep as provi


pdbDIR, pdbPARENT_DIR, pdbTMPL_DIR = _dir_init( __file__, "pdb" )
DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "voronoia" )
VOLUME_CMD = os.path.join( TMPL_DIR, "get_volume_32.exe" )

logging.basicConfig()
LOG = logging.getLogger('voronoia')
LOG.setLevel( logging.ERROR )


HOLE_PARTLY_FILLED = 1
HOLE_NOT_FILLED = 2
HOLE_PARTLY_FILLED_HETS_REMOVED = 3
HOLE_FILLED_HETS_REMOVED = 4

HoleNeighbour = collections.namedtuple( "HoleNeighbour", [
    "atomno", "atomname", "resno", "resname", "chain"
])
VolHole = collections.namedtuple( "VolHole", [
    "no", "type", "neighbours"
])
InfoRecord = collections.namedtuple( "InfoRecord", [
    "pdb_id", "pdb_res", "pdb_title", "pdb_experiment", "pdb_zscorerms"
])

def parse_vol( vol_file, pdb_file ):
    def get_dicts( base ):
        def get_json_dict( tmp_dir, json_tmp_file ):
            with open( os.path.join( tmp_dir, json_tmp_file ), "r" ) as fp2:
                return json.load( fp2 )
        dics = []
        for dic in [
            (pdbTMPL_DIR, "protor_"+base+".json"),
            (TMPL_DIR, "ref_protor_packing_"+base+".json"),
            (TMPL_DIR, "ref_protor_deviation_"+base+".json")
        ]:
            dics.append( get_json_dict(dic[0], dic[1] ))
        return dics
    def dic_lookup( res_sort, residue, atom, atomno, pdb_file ):
        log2_list = ''
        protordic, ref_packing, ref_deviation = get_dicts( res_sort )
        try:
            reffi = float( ref_deviation[residue][protordic[residue][atom]] )
        except Exception:
            reffi = 0.0
        try:
            p_dict[ atomno ] = residue+'|'+ protordic[residue][atom]
        except Exception:
            p_dict[ atomno ] = residue+ '| '
            log2_list += str(pdb_file+'\t'+str(atomno)+ '\t'+residue+ '\t'+atom+'\t protor\n')
        if reffi != 0.0:
            try:
                ref_pa = float( ref_packing[residue][protordic[residue][atom]] )
            except Exception:
                ref_pa = 1.0
            return ( ( packdens-ref_pa )/reffi )**2, log2_list
        else:
            try:
                ref_pa = float( ref_packing[residue][protordic[residue][atom]] )
            except Exception:
                log2_list += str(pdb_file+'\t'+str(atomno)+ '\t'+residue+ '\t'+atom+'\t ref\n')
                ref_pa = 1.0
            return 10, log2_list
    pdb_coord_dict = provi.get_pdb_coord_dict( pdb_file )
    pdb_index_dict = provi.get_pdb_index_dict( pdb_file )
    vol_index_dict = {}

    vol_lines = [None] * len(pdb_coord_dict)
    pd_dict = {}
    buried_dict = {}
    holes = collections.OrderedDict()
    hole_types = []
    hole_list = []
    nrholes = {}
    zs_dict = {}
    p_dict = {}
    log_list = ''
    
    i = 1
    with open( vol_file, "r" ) as fp:
        for l in fp:
            if l.startswith('ATOM') or l.startswith('HETATM'):
                key = ( float(l[30:38]), float(l[38:46]), float(l[46:54]) )
                if key in pdb_coord_dict:
                    j = pdb_coord_dict[key]
                    vol_lines[ j-1 ] = l
                    vol_index_dict[ i ] = j
                    i += 1
                else:
                    raise Exception( 
                        "vol coords not in pdb coords dict. %s %s" % (
                            str(key), l
                        )
                    )
            elif l.startswith('HOLE NUMBER'):
                ls = map( int, l[12:].split() )
                holes[ ls[0] ] = ls[1:]
            elif l.startswith('NRHOLE'):
                p = ( ".* (\d+) \(holes partly filled\)"
                        ".* (\d+) \(holes not filled\)"
                        ".* (\d+) \(partly filled holes with Hets removed\)"
                        ".* (\d+) \(filled holes with Hets removed\).*" )
                m = re.match( p, l )
                if m:
                    nrholes = {
                        "partly_filled": int( m.group(1) ),
                        "not_filled": int( m.group(2) ),
                        "partly_filled_hets_removed": int( m.group(3) ),
                        "filled_hets_removed": int( m.group(4) ),
                    }
                    hole_types = ( 
                        [HOLE_PARTLY_FILLED] * int( m.group(1) ) +
                        [HOLE_NOT_FILLED] * int( m.group(2) ) +
                        [HOLE_PARTLY_FILLED_HETS_REMOVED] * int( m.group(3) ) +
                        [HOLE_FILLED_HETS_REMOVED] * int( m.group(4) )
                    )
                else:
                    raise Exception( "error parsing nrholes record" )

    if len(pdb_coord_dict) != len(vol_index_dict):
        LOG.debug( "number of atoms in vol and pdb coord dicts differs." )

    for no, h in holes.iteritems():
        neighbours = []
        for nb in h:
            if nb in vol_index_dict:
                neighbours.append( vol_index_dict[nb] )
            else:
                raise Exception( "hole neighbour index not found. %s" % nb )
        neighbours.sort()
        neighbours2 = []
        for nb in neighbours:
            pl = pdb_index_dict[ nb ]
            neighbours2.append(
                HoleNeighbour(
                    int( pl[6:11] ),    # atomno
                    pl[12:16],          # atomname
                    int( pl[22:26] ),   # resno
                    pl[17:20],          # resname
                    pl[21:22],          # chain
                )
            )
        hole_list.append( 
            VolHole( no, hole_types[ no-1 ], neighbours2 )
        )
    for i, l in enumerate(vol_lines):
        if l:
            ls = l.split()
            buried = int(ls[-1])    # BURIED
            vdwvol = float(ls[-3])  # VOLUME INSIDE VAN-DER-WAALS SPHERE
            sevol = float(ls[-2])   # VOLUME IN 1.4 ANGSTROM LAYER OUTSIDE VDW-SPHERE
            if vdwvol==0.0:
                packdens = 0.0
                LOG.error( "vdw volume zero. %s" % l )
            elif (vdwvol+sevol)==0.0:
                packdens = 0.0
                LOG.error( 
                    "sum of vdw volume and excluded volume zero. %s" % l 
                )
            else:
                packdens = (vdwvol/(vdwvol+sevol))
            atomno = int( pdb_index_dict[ i+1 ][6:11] )
            pd_dict[ atomno ] = packdens
            buried_dict[ atomno ] = buried
            #calculate z-score
            residue = pdb_index_dict[ i+1 ][17:20].split()[0]
            atom = pdb_index_dict[ i+1 ][12:16].split()[0]
            atom = atom.replace( "\'", "*" )
            zscore_per_atom_rna = 10; zscore_per_atom_prot = 10
            tmp_list=''
            if residue in ['G', 'C', 'A', 'U']:
                zscore_per_atom_rna, tmp_list = dic_lookup( 'rna', residue, atom, atomno, pdb_file )
            else:
                zscore_per_atom_prot, tmp_list = dic_lookup( 'nuc', residue, atom, atomno, pdb_file )
            if zscore_per_atom_rna != 10:
                
                zscore = zscore_per_atom_rna
            else: zscore = zscore_per_atom_prot
            zs_dict[ atomno ] = zscore
            log_list += tmp_list
            
            
        else:
            LOG.debug( 
                "no volume data for line: %s" % (
                    pdb_index_dict[ i+1 ].strip('\n') 
                )
            )
    return {
        "nrholes": nrholes,
        "holes": hole_list,
        "packdens": pd_dict,
        "buried": buried_dict,
        "zscore": zs_dict,
        "protors": p_dict,
        "log_list": log_list
    }

def make_ref( tool_results ):
    def dic_in_dic( elem, ref_dic, out ):
        short = {}
        short[elem.split('|')[1]] = ref_dic
        out[elem.split('|')[0]].update(short)
    ref_dic_dens = collections.defaultdict( list )
    out_dev = collections.defaultdict( dict )
    out_dens = collections.defaultdict( dict )
    log_list = ''
    for t in tool_results:
        log_list += t.log_list
        for elem in t.protor:
            if t.burried[elem] and not (t.protor[elem].split('|')[0] in ['G', 'C', 'A', 'U']):
                ref_dic_dens[ t.protor[elem] ].append( t.packdens[elem] )
    for elem in ref_dic_dens:
        dic_in_dic( elem, std(ref_dic_dens[elem]), out_dev )
        ref_dic_dens[elem]=sum(ref_dic_dens[elem])/len(ref_dic_dens[elem])
        dic_in_dic( elem, ref_dic_dens[elem], out_dens )
    return out_dens, out_dev, log_list


def make_nrhole_pdb( pdb_input, holes, nh_file, std_file, pymol_file):
    npdb = numpdb.NumPdb( pdb_input )
    sele2={'record':'HETATM'}
    last_hetatm=int(npdb.get('atomno', **sele2)[-1])
    fi=open(nh_file, 'w')
    preline=""
    std_dct={}
    neighbours={}
    with open(pdb_input, 'r') as fp:
        for line in fp:
            if preline!=line[0:6] and preline=="HETATM":
                for hno, hole in enumerate(holes):
                    chain=''
                    xyz_list=[]
                    atom_list=[]
                    for hn in hole[2]:
                        atomno, atomname, resno, resname, chain=hn
                        sele={'atomno':atomno, 'chain':chain}
                        xyz=npdb.get('xyz', **sele)[0]
                        xyz_list.append(xyz)
                        atom_list.append([atomno, resno, chain])
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
                    std_dct[last_hetatm+hno]=mean
                    fi.write(numpdb.pdb_line( natom ))
                fi.write(line)
            else:
                fi.write(line)
            preline=line[0:6]
    with open( std_file, "w" ) as fp:
        json.dump( std_dct, fp )
    with open(pymol_file, 'w') as fp:
        code='import pymol; \n'+\
            'pymol.finish_launching() \n'+\
            'neighbours='+json.dumps(neighbours)+' \n'+\
            'mean_dic='+json.dumps(std_dct)+' \n'+\
            'pymol.cmd.load("'+nh_file+'") \n'+\
            'pymol.cmd.hide( representation="line") \n'+\
            'pymol.cmd.show( representation="cartoon") \n'+\
            'pymol.cmd.select( "neh", "resname NEH") \n'+\
            'pymol.cmd.show(representation="spheres", selection="neh") \n'+\
            'pymol.cmd.spectrum( "b", "blue_white_red", "neh") \n'+\
            'for index, elem in enumerate(mean_dic): \n'+\
            '    pymol.cmd.select( "temp", "(resi "+str(index)+" and resname NEH)" ) \n'+\
            '    pymol.cmd.alter("temp", "vdw="+str(mean_dic[elem]/2)) \n'+\
            '    pymol.cmd.rebuild() \n'+\
            'for index, elem in enumerate(neighbours): \n'+\
            '    liste="" \n'+\
            '    for neighbour in neighbours[elem]: \n'+\
            '        if liste!="": \n'+\
            '            liste=liste+"+" \n'+\
            '        liste=liste+"(id "+str(neighbour[0])+" and resi "+str(neighbour[1])+" and chain "+neighbour[2]+")" \n'+\
            '    pymol.cmd.select( "neighbour_"+str(elem), liste ) \n'
        fp.write(code)


# get_volume.exe ex:0.1 rad:protor i:file.pdb o:out.vol
class Voronoia( CmdTool, ProviMixin, ParallelMixin, RecordsMixin ):
    """A wrapper around the 'voronoia' aka 'get_volume' programm."""
    args = [
        _( "pdb_input", type="file", ext="pdb" ),
        _( "ex", type="float", range=[0.01, 0.5], step=0.01, default=0.1 ),
        _( "radii", type="str", options=["protor"], default="protor" ),
        # TODO run the tool multiple times 
        # when it fails without the shuffle option
        _( "shuffle", type="bool", default=False, 
            help="slightly changes the input coordinates to"
                "circumvent numerical problems" ),
        _( "make_reference|mr", type="bool", default=False ),
        _( "analyze_only|ao", type="bool", default=False ),
        _( "get_nrholes|gh", type="bool", default=False )
    ]
    out = [
        _( "vol_file", file="{pdb_input.stem}.vol" ),
        _( "log_file", file="{pdb_input.stem}.log" ),
        _( "dens_file", file="ref_protor_packing.json", optional=True ),
        _( "dev_file", file="ref_protor_deviation.json", optional=True ),
        _( "protor_log_file", file="protor.log", optional=True ),
        _( "nh_file", file="{pdb_input.stem}_nh.pdb", optional=True ),
        _( "pymol_file", file="pymol_settings.py", optional=True ),
        _( "std_file", file="{pdb_input.stem}_std.json", optional=True ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "voronoia.provi"
    RecordsClass = InfoRecord
    def _init( self, *args, **kwargs ):
        self._init_records( None, **kwargs )
        self._init_parallel( self.pdb_input, **kwargs )
        if not self.analyze_only:
            if not self.parallel: #and not self.analyze_only:
                self.cmd = [ 
                    "wine", VOLUME_CMD, 
                    "ex:%0.1f"%float(self.ex), 
                    "rad:%s"%self.radii,
                    "x:yes",
                    "l:%s" %self.log_file,
                    "i:%s"%self.pdb_input, 
                    "o:%s"%self.vol_file
                ]
                if self.shuffle:
                    self.cmd.append( "sh:y" )
            else:
                self.cmd = None
        else:
            self.cmd = None
    def _pre_exec( self ):
        utils.path.remove( self.vol_file )
        utils.path.remove( self.log_file )
    def _post_exec( self ):
        if not self.parallel:
            provi.prep_volume( self.vol_file, self.pdb_input )
            self._make_provi_file(
                pdb_file=self.relpath( self.pdb_input ),
                vol_file=self.relpath( self.vol_file )
            )
            self.info = numpdb.NumPdb( self.pdb_input, features={
                "phi_psi": False, 
                "info": True,
                "backbone_only": True
            })._info
            #get info for InfoRecord
            dicts = self.get_vol()
            self.zscores = dicts['zscore']
            self.protor = dicts['protors']
            self.burried = dicts['buried']
            self.packdens = dicts['packdens']
            self.log_list = dicts['log_list']
            self.holes = dicts['holes']
            zscore_allb = []; zscore_all = 0
            for elem in self.zscores:
                if self.burried[elem]:
                    zscore_allb.append(self.zscores[elem])
            zscore_all = sum(zscore_allb)
            pdbid, ext = os.path.splitext(os.path.basename(self.pdb_input))
            self.zscorerms = zscore_all/len(zscore_allb)
            self.records = [
                InfoRecord(
                    pdbid, self.info["resolution"],
                    self.info["title"], self.info["experiment"],
                    self.zscorerms
                )
            ]
            self.write()
            if self.get_nrholes:
                make_nrhole_pdb(self.pdb_input,self.holes, self.nh_file, self.std_file, self.pymol_file)
        if self.parallel and self.make_reference:
            dict_dens, dict_dev, log_list = make_ref( self.tool_results )
            d = ( self.dens_file, dict_dens ), ( self.dev_file, dict_dev )
            for fname, dct in d:
                with open( fname, "w" ) as fp:
                    json.dump( dct, fp )
            with open(self.protor_log_file, 'w') as fp:
                fp.write( log_list )
    @memoize_m
    def get_vol( self ):
        return parse_vol( self.vol_file, self.pdb_input )






