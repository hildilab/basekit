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
from utils.tool import RecordsMixin, PyTool, SqliteBackend
from utils.listing import merge_dic_list
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
VoronoiaDbRecord = collections.namedtuple( "VoronoiaDbRecord", [
    "pdb_id", "pdb_title", "pdb_keywords", "pdb_experiment",
   # "pdb_header",
    "pdb_resolution",
    "pdb_zscorerms",
    "curated_representative", "curated_related",
    "status",
    "tm_packdens_protein_buried", "tm_water_count", "tm_residue_count", 
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
    pd_at_dict = {} #(res, atom):[pd, pd,...]
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
            residue = pdb_index_dict[ i+1 ][17:20].split()[0]
            atom = pdb_index_dict[ i+1 ][12:16].split()[0]
            atom = atom.replace( "\'", "*" )
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
                if residue in pd_at_dict:
                    if atom in pd_at_dict[residue]:
                        pd_at_dict[residue][atom]=pd_at_dict[residue][atom]+[packdens]
                    else:
                        pd_at_dict[residue][atom]=[packdens]
                else:
                    pd_at_dict[residue]={}
                    pd_at_dict[residue][atom]=[packdens]
            atomno = int( pdb_index_dict[ i+1 ][6:11] )
            pd_dict[ atomno ] = packdens
            buried_dict[ atomno ] = buried
            #calculate z-score
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
        "log_list": log_list,
        "pd_at_dict": pd_at_dict
    }

def make_ref( tool_results ):
    def dic_in_dic( elem, ref_dic, out ):
        short = {}
        short[elem.split('|')[1]] = ref_dic
        out[elem.split('|')[0]].update(short)
    
    ref_dic_dens = collections.defaultdict( list )
    out_dev = collections.defaultdict( dict )
    out_dens = collections.defaultdict( dict )
    out_pd_at_dict = {}
    log_list = ''
    for t in tool_results:
        log_list += t.log_list
        for elem in t.protor:
            if t.burried[elem] and not (t.protor[elem].split('|')[0] in ['G', 'C', 'A', 'U']):
                ref_dic_dens[ t.protor[elem] ].append( t.packdens[elem] )
        out_pd_at_dict = merge_dic_list(out_pd_at_dict, t.pd_at_dict)
    for elem in ref_dic_dens:
        dic_in_dic( elem, std(ref_dic_dens[elem]), out_dev )
        ref_dic_dens[elem]=sum(ref_dic_dens[elem])/len(ref_dic_dens[elem])
        dic_in_dic( elem, ref_dic_dens[elem], out_dens )
    return out_dens, out_dev, log_list, out_pd_at_dict

def holes_mean_std_helper(hno, last_hetatm, chain, xyz_list, resno):
    
    xyz2=xyz_list
    coords = np.mean(xyz_list, axis=0)
    xyz_list_dist=[]
    for co in xyz2:
        xyz_list_dist.append(np.sqrt(np.sum((coords-co)**2)))
    mean_num=np.mean(xyz_list_dist, axis=0)
    std_num=np.std(xyz_list_dist)
    natom={
        "record":"HETATM","atomno": last_hetatm+hno, "atomname": " CA ",
        "resname": "NEH", "chain": chain, "resno": resno,
        "x": coords[0], "y": coords[1], "z": coords[2], "bfac": std_num, "element": " C"
    }
    return numpdb.pdb_line( natom ), mean_num

def make_nrhole_pdb( pdb_input, holes, nh_file, mean_file):
    features={"protein_only": False}
    npdb = numpdb.NumPdb( pdb_input, features=features )
    try:
        sele2={'record':'HETATM'}
        last_hetatm=int(npdb.get('atomno', **sele2)[-1])
        last_hetresno=int(npdb.get('resno', **sele2)[-1])  
    except:
        sele3={'record':'ATOM  '}
        last_hetatm=int(npdb.get('atomno', **sele3)[-1])
        last_hetresno=int(npdb.get('resno', **sele3)[-1])
    fi=open(nh_file, 'w')
    preline=""
    mean_dct={}; neighbours={}; mean_lst=[]
    exiting=False; writing=False
    with open(pdb_input, 'r') as fp:
        for line in fp:
            if str(last_hetatm) in line[0:13]:
                writing=True
                fi.write(line)
            elif writing==True and exiting==False:
                for hno, hole in enumerate(holes):
                    chain=''
                    xyz_list=[]
                    atom_list=[]
                    for hn in hole[2]:
                        try:
                            atomno, atomname, resno, resname, chain=hn
                            sele={'atomno':atomno, 'chain':chain}
                            xyz=npdb.get('xyz', **sele)[0]
                            xyz_list.append(xyz)
                            atom_list.append([atomno, resno, chain])
                        except:
                            break
                    neighbours[hno]=atom_list
                    new_line, mean_num = holes_mean_std_helper(hno, last_hetatm, chain, xyz_list, resno)
                    fi.write(new_line)
                    mean_dct[last_hetatm+hno]=mean_num
                    mean_lst.append(str(atomno)+'_'+str(resno)+'_'+chain)
                fi.write(line)
                exiting=True
            else:
                fi.write(line)
            preline=line[0:6]
    with open( mean_file, "w" ) as fp:
        json.dump( mean_dct, fp )
    return neighbours, mean_dct, last_hetresno, mean_lst


###-->from old server

#packing density + extended_vol_file
def packing_density_file(volfile):
    fil= volfile
    infile = open(fil,"r")
    outfile = file(fil+'.extended_temp.vol', "w")
    z_score_all=0
    z_score_all_protein=0
    z_score_all_atNum=0
    z_score_all_atNum_protein=0
    for line in infile:
        if (line.startswith('ATOM') or line.startswith('HETATM')):
            m=line.split()
            vdwvol=float(m[-3])  #VOLUME INSIDE VAN-DER-WAALS SPHERE
            sevol=float(m[-2])   #VOLUME IN 1.4 ANGSTROM LAYER OUTSIDE VDW-SPHERE
            if vdwvol==0.0:
                ppdd='0.0'
               # print fil, line
            elif (vdwvol+sevol)==0.0:
                ppdd= '999.99'
                #print fil, line
            else:
                ppdd="%.3f" % (vdwvol/(vdwvol+sevol))
            if ((line.startswith('ATOM') and (line.endswith('1\r\n') or line.endswith('1\n'))) or (line.startswith('HETATM') and (line.endswith('1\r\n') or line.endswith('1\n')))):
                atom=line[12:16].split()[0]
                residue=line[17:20].split()[0]
                m=line.split()
                vdwvol=m[-3]  #VOLUME INSIDE VAN-DER-WAALS SPHERE
                sevol=m[-2]   #VOLUME IN 1.4 ANGSTROM LAYER OUTSIDE VDW-SPHERE
               # print vdwvol, sevol
                if float(vdwvol)==0.0:
                    zscore_per_atom_str=""
                elif (float(vdwvol)+float(sevol))==0.0:
                    zscore_per_atom_str=""
                else:
                    density=(float(vdwvol)/(float(vdwvol)+float(sevol)))
                    if residue in ['G', 'C', 'A', 'U']:
                        protordic={"G":{"C1*":"C4H1b","C3*":"C4H1s","O6":"O1H0u","C5*":"C4H2s","O4*":"O2H0s","O1P":"O1H0u","C8":"C3H1s","O2*":"O2H1u","C2":"C3H0s","C6":"C3H0s","C5":"C3H0s","C4":"C3H0s","O2P":"O1H0u","P":"P4H0u","C2*":"C4H1s","N1":"N3H1s","N2":"N3H2u","N3":"N2H0b","C4*":"C4H1b","N7":"N2H0b","O5*":"O2H0b","N9":"N3H0u","O3*":"O2H0b"},"A":{"C1*":"C4H1b","C3*":"C4H1s","C5*":"C4H2s","O4*":"O2H0s","O1P":"O1H0u","C8":"C3H1s","O2*":"O2H1u","N6":"N3H2u","C2":"C3H1s","C6":"C3H0s","C5":"C3H0s","C4":"C3H0s","O2P":"O1H0u","P":"P4H0u","C2*":"C4H1s","N1":"N2H0s","N3":"N2H0b","C4*":"C4H1b","N7":"N2H0b","O5*":"O2H0b","N9":"N3H0u","O3*":"O2H0b"},"C":{"C5*":"C4H2s","C1*":"C4H1b","O4*":"O2H0s","C3*":"C4H1s","O5*":"O2H0b","O1P":"O1H0u","O2P":"O1H0u","C6":"C3H1t","P":"P4H0u","O3*":"O2H0b","O2*":"O2H1u","C4*":"C4H1b","C2":"C3H0s","C2*":"C4H1s","N1":"N3H0u","N3":"N2H0s","N4":"N3H2u","O2":"O1H0u","C5":"C3H1b","C4":"C3H0s"},"U":{"C5*":"C4H2s","C1*":"C4H1b","O4*":"O2H0s","N3":"N3H1s","C3*":"C4H1s","O5*":"O2H0b","O1P":"O1H0u","O2P":"O1H0u","C6":"C3H1t","P":"P4H0u","O3*":"O2H0b","O2*":"O2H1u","C4*":"C4H1b","C2":"C3H0s","C2*":"C4H1s","N1":"N3H0u","O4":"O1H0u","O2":"O1H0u","C5":"C3H1b","C4":"C3H0s"},"T":{"C5*":"C4H2s","C1*":"C4H1b","O4*":"O2H0s","N3":"N3H1s","C3*":"C4H1s","O5*":"O2H0b","O1P":"O1H0u","O2P":"O1H0u","C6":"C3H1t","P":"P4H0u","O3*":"O2H0b","O2*":"O2H1u","C4*":"C4H1b","C2":"C3H0s","C2*":"C4H1s","N1":"N3H0u","O4":"O1H0u","C7":"C4H3u","O2":"O1H0u","C5":"C3H0b","C4":"C3H0s"}}
                        ref_packing={"G":{"C4H1b":"0.748","C4H1s":"0.805","O1H0u":"0.567","C4H2s":"0.673","O2H0s":"0.699","C3H1s":"0.672","O2H1u":"0.572","C3H0s":"0.829","P4H0u":"0.807","N3H1s":"0.765","N3H2u":"0.620","N2H0b":"0.664","O2H0b":"0.702","N3H0u":"0.854"},"A":{"C4H1b":"0.754","C4H1s":"0.805","C4H2s":"0.669","O2H0s":"0.689","O1H0u":"0.559","C3H1s":"0.673","O2H1u":"0.569","N3H2u":"0.613","C3H0s":"0.833","P4H0u":"0.808","N2H0s":"0.744","N2H0b":"0.662","O2H0b":"0.703","N3H0u":"0.853"},"C":{"C4H2s":"0.667","C4H1b":"0.740","O2H0s":"0.693","C4H1s":"0.808","O2H0b":"0.701","O1H0u":"0.567","C3H1t":"0.707","P4H0u":"0.806","O2H1u":"0.573","C3H0s":"0.831","N3H0u":"0.849","N2H0s":"0.762","N3H2u":"0.603","C3H1b":"0.635"},"U":{"C4H2s":"0.676","C4H1b":"0.750","O2H0s":"0.689","N3H1s":"0.736","C4H1s":"0.801","O2H0b":"0.700","O1H0u":"0.550","C3H1t":"0.708","P4H0u":"0.813","O2H1u":"0.579","C3H0s":"0.809","N3H0u":"0.845","C3H1b":"0.631"}}
                        ref_deviation={"G":{"C4H1b":"0.090","C4H1s":"0.087","O1H0u":"0.083","C4H2s":"0.087","O2H0s":"0.089","C3H1s":"0.091","O2H1u":"0.089","C3H0s":"0.092","P4H0u":"0.089","N3H1s":"0.084","N3H2u":"0.084","N2H0b":"0.087","O2H0b":"0.088","N3H0u":"0.085"},"A":{"C4H1b":"0.074","C4H1s":"0.065","C4H2s":"0.064","O2H0s":"0.064","O1H0u":"0.064","C3H1s":"0.076","O2H1u":"0.086","N3H2u":"0.060","C3H0s":"0.000","P4H0u":"0.086","N2H0s":"0.064","N2H0b":"0.064","O2H0b":"0.063","N3H0u":"0.060"},"C":{"C4H2s":"0.073","C4H1b":"0.077","O2H0s":"0.083","C4H1s":"0.073","O2H0b":"0.083","O1H0u":"0.069","C3H1t":"0.078","P4H0u":"0.083","O2H1u":"0.083","C3H0s":"0.078","N3H0u":"0.068","N2H0s":"0.070","N3H2u":"0.068","C3H1b":"0.079"},"U":{"C4H2s":"0.087","C4H1b":"0.088","O2H0s":"0.089","N3H1s":"0.08","C4H1s":"0.087","O2H0b":"0.089","O1H0u":"0.087","C3H1t":"0.089","P4H0u":"0.089","O2H1u":"0.090","C3H0s":"0.088","N3H0u":"0.086","C3H1b":"0.089"}}
                        if atom=='OP1':atom='O1P'
                        if atom=='OP2':atom='O2P'
                        if atom=='OP3':atom='O2P'
                        for elem in atom:
                            if elem=='\'': atom=atom.replace("\'", "*")
                        reffi=float(ref_deviation[residue][protordic[residue][atom]])
                        if reffi==0.0:
                            zscore_per_atom_str=""
                        else:
                            zscore_per_atom=((density-float(ref_packing[residue][protordic[residue][atom]]))/float(ref_deviation[residue][protordic[residue][atom]]))**2
                            z_score_all=z_score_all+zscore_per_atom
                            zscore_per_atom_str="%.3f" % zscore_per_atom
                            z_score_all_atNum=z_score_all_atNum+1
                    elif residue in ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]:
                        ref_packing_protein={"ALA":{"C":"0.810","CA":"0.756","CB":"0.582","N":"0.754","O":"0.606"},"ARG":{"C":"0.820","CA":"0.775","CB":"0.664","CD":"0.674","CG":"0.660","CZ":"0.778","N":"0.763","NE":"0.696","NH1":"0.616","NH2":"0.601","O":"0.610"},"ASN":{"C":"0.813","CA":"0.779","CB":"0.666","CG":"0.777","N":"0.754","ND2":"0.604","O":"0.612","OD1":"0.587"},"ASP":{"C":"0.816","CA":"0.775","CB":"0.656","CG":"0.742","N":"0.747","O":"0.614","OD1":"0.616","OD2":"0.602"},"CYS":{"C":"0.812","CA":"0.774","CB":"0.661","N":"0.747","O":"0.601","SG":"0.592"},"GLN":{"C":"0.820","CA":"0.776","CB":"0.663","CD":"0.774","CG":"0.664","N":"0.769","NE2":"0.602","O":"0.616","OE1":"0.569"},"GLU":{"C":"0.822","CA":"0.775","CB":"0.659","CD":"0.738","CG":"0.653","N":"0.766","O":"0.612","OE1":"0.597","OE2":"0.590"},"GLY":{"C":"0.781","CA":"0.660","N":"0.728","O":"0.598"},"HIS":{"C":"0.819","CA":"0.779","CB":"0.660","CD2":"0.633","CE1":"0.629","CG":"0.849","N":"0.758","ND1":"0.681","NE2":"0.677","O":"0.606"},"ILE":{"C":"0.842","CA":"0.791","CB":"0.754","CD1":"0.571","CG1":"0.651","CG2":"0.588","N":"0.762","O":"0.611"},"LEU":{"C":"0.819","CA":"0.786","CB":"0.659","CD1":"0.572","CD2":"0.571","CG":"0.738","N":"0.763","O":"0.613"},"LYS":{"C":"0.819","CA":"0.777","CB":"0.658","CD":"0.645","CE":"0.653","CG":"0.665","N":"0.761","NZ":"0.643","O":"0.612"},"MET":{"C":"0.819","CA":"0.778","CB":"0.662","CE":"0.597","CG":"0.664","N":"0.765","O":"0.609","SD":"0.606"},"PHE":{"C":"0.823","CA":"0.777","CB":"0.657","CD1":"0.634","CD2":"0.624","CE1":"0.606","CE2":"0.603","CG":"0.861","CZ":"0.605","N":"0.759","O":"0.605"},"PRO":{"C":"0.815","CA":"0.767","CB":"0.632","CD":"0.678","CG":"0.632","N":"0.853","O":"0.603"},"SER":{"C":"0.809","CA":"0.768","CB":"0.642","N":"0.739","O":"0.608","OG":"0.595"},"THR":{"C":"0.824","CA":"0.782","CB":"0.730","CG2":"0.584","N":"0.752","O":"0.611","OG1":"0.596"},"TRP":{"C":"0.828","CA":"0.775","CB":"0.658","CD1":"0.627","CD2":"0.788","CE2":"0.807","CE3":"0.638","CG":"0.839","CH2":"0.612","CZ2":"0.618","CZ3":"0.608","N":"0.758","NE1":"0.649","O":"0.611"},"TYR":{"C":"0.824","CA":"0.779","CB":"0.657","CD1":"0.637","CD2":"0.631","CE1":"0.620","CE2":"0.619","CG":"0.860","CZ":"0.787","N":"0.761","O":"0.606","OH":"0.572"},"VAL":{"C":"0.835","CA":"0.787","CB":"0.735","CG1":"0.585","CG2":"0.586","N":"0.757","O":"0.609"}}
                        ref_deviation_protein={"ALA":{"C":"0.046","CA":"0.053","CB":"0.049","N":"0.058","O":"0.058"},"ARG":{"C":"0.050","CA":"0.051","CB":"0.049","CD":"0.055","CG":"0.053","CZ":"0.057","N":"0.053","NE":"0.063","NH1":"0.062","NH2":"0.065","O":"0.060"},"ASN":{"C":"0.050","CA":"0.050","CB":"0.051","CG":"0.053","N":"0.055","ND2":"0.066","O":"0.059","OD1":"0.074"},"ASP":{"C":"0.047","CA":"0.051","CB":"0.055","CG":"0.055","N":"0.056","O":"0.061","OD1":"0.075","OD2":"0.075"},"CYS":{"C":"0.049","CA":"0.052","CB":"0.051","N":"0.055","O":"0.058","SG":"0.067"},"GLN":{"C":"0.049","CA":"0.049","CB":"0.052","CD":"0.053","CG":"0.051","N":"0.052","NE2":"0.058","O":"0.062","OE1":"0.072"},"GLU":{"C":"0.047","CA":"0.049","CB":"0.051","CD":"0.057","CG":"0.053","N":"0.052","O":"0.058","OE1":"0.074","OE2":"0.080"},"GLY":{"C":"0.056","CA":"0.052","N":"0.058","O":"0.060"},"HIS":{"C":"0.052","CA":"0.050","CB":"0.053","CD2":"0.062","CE1":"0.062","CG":"0.044","N":"0.056","ND1":"0.077","NE2":"0.090","O":"0.058"},"ILE":{"C":"0.049","CA":"0.051","CB":"0.052","CD1":"0.053","CG1":"0.055","CG2":"0.049","N":"0.048","O":"0.057"},"LEU":{"C":"0.049","CA":"0.045","CB":"0.048","CD1":"0.050","CD2":"0.051","CG":"0.058","N":"0.052","O":"0.060"},"LYS":{"C":"0.049","CA":"0.052","CB":"0.051","CD":"0.055","CE":"0.060","CG":"0.051","N":"0.053","NZ":"0.079","O":"0.060"},"MET":{"C":"0.050","CA":"0.050","CB":"0.053","CE":"0.053","CG":"0.053","N":"0.053","O":"0.060","SD":"0.058"},"PHE":{"C":"0.058","CA":"0.048","CB":"0.051","CD1":"0.054","CD2":"0.055","CE1":"0.059","CE2":"0.058","CG":"0.039","CZ":"0.060","N":"0.055","O":"0.057"},"PRO":{"C":"0.042","CA":"0.055","CB":"0.055","CD":"0.051","CG":"0.055","N":"0.043","O":"0.064"},"SER":{"C":"0.048","CA":"0.053","CB":"0.053","N":"0.055","O":"0.061","OG":"0.069"},"THR":{"C":"0.051","CA":"0.052","CB":"0.054","CG2":"0.057","N":"0.052","O":"0.061","OG1":"0.075"},"TRP":{"C":"0.055","CA":"0.047","CB":"0.047","CD1":"0.054","CD2":"0.051","CE2":"0.057","CE3":"0.057","CG":"0.043","CH2":"0.058","CZ2":"0.054","CZ3":"0.058","N":"0.056","NE1":"0.060","O":"0.061"},"TYR":{"C":"0.058","CA":"0.048","CB":"0.050","CD1":"0.054","CD2":"0.054","CE1":"0.056","CE2":"0.057","CG":"0.040","CZ":"0.056","N":"0.054","O":"0.058","OH":"0.058"},"VAL":{"C":"0.050","CA":"0.050","CB":"0.057","CG1":"0.052","CG2":"0.051","N":"0.048","O":"0.057"}}
                        if ref_packing_protein[residue].has_key(atom):
                            for elem in atom:
                                if elem=='\'': atom=atom.replace("\'", "*")
                            reffi=float(ref_deviation_protein[residue][atom])
                            if reffi==0.0:
                                zscore_per_atom_str=""
                            else:
                                zscore_per_atom_protein=((density-float(ref_packing_protein[residue][atom]))/float(ref_deviation_protein[residue][atom]))**2
                                z_score_all_protein=z_score_all_protein+zscore_per_atom_protein
                                zscore_per_atom_str="%.3f" % zscore_per_atom_protein
                                z_score_all_atNum_protein=z_score_all_atNum_protein+1
                        else: zscore_per_atom_str=""
                    else: zscore_per_atom_str=""
            else: zscore_per_atom_str=""
            newline=line[:82]+' '+ppdd+' '+zscore_per_atom_str+'\r\n'
            outfile.write(newline)
        else:
            outfile.write(line)
    infile.close()
    outfile.close()
    infilesec=open(fil+'.extended_temp.vol',"r")
    outfilesec=file(fil+'.extended.vol', "w")
    
    for line in infilesec:
        if line.startswith('RESOLT') :
            outfilesec.write(line)
            if z_score_all_atNum!=0:
                outfilesec.write('ZSCORE_RMS_RNA '+"%.3f" % (math.sqrt(z_score_all/z_score_all_atNum))+'\r\n')
                if z_score_all_atNum_protein!=0:
                    outfilesec.write('ZSCORE_RMS_ALL '+"%.3f" % (math.sqrt(z_score_all_protein+z_score_all/z_score_all_atNum_protein+z_score_all_atNum))+'\r\n')
                else: outfilesec.write('ZSCORE_RMS_ALL '+"%.3f" % (math.sqrt(z_score_all/z_score_all_atNum))+'\r\n')
            else: outfilesec.write('ZSCORE_RMS '+"%.3f" % (math.sqrt(z_score_all/0.00001))+'\r\n')
            if z_score_all_atNum_protein!=0:
                outfilesec.write('ZSCORE_RMS_PROTEIN '+"%.3f" % (math.sqrt(z_score_all_protein/z_score_all_atNum_protein))+'\r\n')
                if z_score_all_atNum!=0:
                    outfilesec.write('ZSCORE_RMS_ALL '+"%.3f" % (math.sqrt(z_score_all_protein+z_score_all/z_score_all_atNum_protein+z_score_all_atNum))+'\r\n')
                else: outfilesec.write('ZSCORE_RMS_ALL '+"%.3f" % (math.sqrt(z_score_all_protein/z_score_all_atNum_protein))+'\r\n')
            else: outfilesec.write('ZSCORE_RMS_PROTEIN '+"%.3f" % (math.sqrt(z_score_all_protein/0.00001))+'\r\n')
        else: outfilesec.write(line)
    infilesec.close()
    os.remove(fil+'.extended_temp.vol')

###from old server <--




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
        _( "get_nrholes|gh", type="bool", default=False ),
        _( "prefix|pre", type="str", default="" )
    ]
    out = [
        _( "vol_file", file="{pdb_input.stem}.vol" ),
        _( "log_file", file="{pdb_input.stem}.log" ),
        _( "dens_file", file="ref_protor_packing.json", optional=True ),
        _( "dev_file", file="ref_protor_deviation.json", optional=True ),
        _( "pd_at_file", file="pd_at_res.json", optional=True ),
        _( "protor_log_file", file="protor.log", optional=True ),
        _( "nh_file", file="{pdb_input.stem}_nh.pdb", optional=True ),
        _( "pymol_file", file="pymol_settings.py", optional=True ),
        _( "mean_file", file="{pdb_input.stem}_mean.json", optional=True ),
        _( "voro_info", file="voronoia_records.json" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "voronoia.provi"
    pymol_tmpl = "pymol_settings.py"
    RecordsClass = VoronoiaDbRecord
    def _init( self, *args, **kwargs ):
        #if self.pdb_input.endswith("pdb"):
        #    self.log("no pdb")
        #    return
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
            
            #get info for InfoRecord
            dicts = self.get_vol()
            self.zscores = dicts['zscore']
            self.protor = dicts['protors']
            self.burried = dicts['buried']
            self.packdens = dicts['packdens']
            self.log_list = dicts['log_list']
            self.holes = dicts['holes']
            self.pd_at_dict = dicts['pd_at_dict']
            with open( self.pd_at_file, "w" ) as fp:
                json.dump(self.pd_at_dict , fp, indent=4 )
            zscore_allb = []; zscore_all = 0
            for elem in self.zscores:
                if self.burried[elem]:
                    zscore_allb.append(self.zscores[elem])
            zscore_all = sum(zscore_allb)
            pdbid, ext = os.path.splitext(os.path.basename(self.pdb_input))
            self.zscorerms = zscore_all/len(zscore_allb)
            try:
                self.info = numpdb.NumPdb( self.pdb_input, features={
                    "phi_psi": False, 
                    "info": True,
                    "backbone_only": True
                })._info
                self.records = [
                    VoronoiaDbRecord(
                        pdbid, self.info["title"], "no data", self.info["experiment"],
                        self.info["resolution"],
                        self.zscorerms,
                        0, 0,
                        "included",
                        0, 0, 0,
                    )
                ]
            except:
                self.records = [
                    VoronoiaDbRecord(
                        pdbid,  "no data", "no data", "no data",
                        0.0,
                        self.zscorerms,
                        0, 0,
                        "included",
                        0, 0, 0,
                    )
                ]
            db = SqliteBackend( "voronoia.sqlite", VoronoiaDbRecord )
            db.write( self.records )
            self.write()
            
            # get the nrholes and the pymol script
            if self.get_nrholes:
                neighbours, mean_dct, last_hetresno, mean_lst = make_nrhole_pdb(self.pdb_input,self.holes, self.nh_file, self.mean_file)
                values_dict={
                    'neighbours':neighbours, 'mean_dct':mean_dct,
                    'nh_file':self.nh_file, 'pref':self.prefix,
                    'last_hetresno': last_hetresno, 'mean_list':mean_lst
                }
                self._make_file_from_tmpl(self.pymol_tmpl, **values_dict)
            
            ###-->added from old server
            
            packing_density_file(self.vol_file)
            ###added form old server<--
        if self.parallel and self.make_reference:
            dict_dens, dict_dev, log_list, out_pd_at_dict = make_ref( self.tool_results )
            d = ( self.dens_file, dict_dens ), ( self.dev_file, dict_dev ), (self.pd_at_file, out_pd_at_dict)
            for fname, dct in d:
                with open( fname, "w" ) as fp:
                    json.dump( dct, fp, indent=4 )
            with open(self.protor_log_file, 'w') as fp:
                fp.write( log_list )
            db = SqliteBackend( "voronoia.sqlite", VoronoiaDbRecord )
            db.write( self.records )
    @memoize_m
    def get_vol( self ):
        return parse_vol( self.vol_file, self.pdb_input )

