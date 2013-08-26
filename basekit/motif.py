from __future__ import with_statement
from __future__ import division




import numpy as np
from numpy import array
from utils.numpdb import NumPdb, numdist
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.pyplot import figure

import string, os, sys, argparse, itertools

from time import gmtime, strftime
import functools, operator, collections
import basekit.utils.numpdb as numpdb
from utils import bio
from collections import defaultdict

import utils.math, utils.path
from utils import (
    try_int, get_index, boolean, iter_window, memoize, 
    copy_dict, dir_walker
)
from pdb import pdb_download
from utils.timer import Timer
from utils.db import get_pdb_files, create_table
from utils.tool import (
    _, PyTool, DbTool, RecordsMixin, ParallelMixin, ProviMixin
)
import utils.path
from cStringIO import StringIO
import utils.numpdb as numpdb
import json
from basekit.msa import Muscle
import linecache
import scipy





DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "motif" )



MotifRecord = collections.namedtuple( 'MotifRecord', [
    'motif', 'pdb_id', 'chain', 'resno', 'rmsd',
    'phi_angles', 'psi_angles', 'l2f_str', 'm2f_str',  'no'
])

StrucRecord = collections.namedtuple( 'StrucRecord', [
    'pdb_id', 'chain', 'resno', 'phi_angles',
    'psi_angles', 'distances', 'no'
])

SuperposeRecord = collections.namedtuple( 'SuperposeRecord', [ ])


class MotifTest(object):
    length = 0 # to be set by classes
    def get_dihedral(self, name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    def make_list(self, *numatoms):
        Ccap = numatoms[3]
        return (
            Ccap['resno'][0],Ccap['chain'][0],(
                self.get_dihedral( "phi", *numatoms ),
                self.get_dihedral( "psi", *numatoms )
            )
        )
    def test(self, *atoms):
        # to be implemented be sub classes
        pass
    

class AlphaBetaMotifTest(MotifTest):
    """
    * H-bond between CO of residue(i)=C3 and NH of residue(i+4)=Cd
    * H-bond between CO of residue(i)=C3 and NH of residue(i+3)=Ccap
    * phi angles of residues(i+1)=C2, (i+2)=C1, (i+3)=Ccap are negative
    Reference:
        Adel Golovin, Kim Henrick;
            MSDmotif: exploring protein sites and motifs (2008), BMC Bioinformatics 9:312
            Appendix A. Small 3D structural motifs
    """
    length = 6
    def test(self, C3, C2, C1, Ccap, Cd, Cdd):
        C3_dist = C3.copy(atomname="O")
        Ccap_dist = Ccap.copy(atomname="N")
        Cd_dist = Cd.copy(atomname="N")
        distance1 = numpdb.numdist( C3_dist, Cd_dist )
        distance2 = numpdb.numdist( C3_dist, Ccap_dist )
        return distance1 < 3.9 and distance2 < 3.9 and C2['phi'][0]<0 and C1['phi'][0]<0 and Ccap['phi'][0]<0

    
class AsxMotifTest():
    """
    * Residue(i)=C3 is Aspartate(D) or Asparagine(Asx)
    * H-bond between (side-chain) O of residue(i)=C3 and (main-chain) NH of residue(i+2)=C1 or residue(i+3)=Ccap
    * H-bond between CO of residue(i)=C3 and NH of residue(i+3)=Ccap or residue(i+4)=Cd
    Reference:
        Adel Golovin, Kim Henrick;
            MSDmotif: exploring protein sites and motifs (2008), BMC Bioinformatics 9:312
            Appendix A. Small 3D structural motifs
    """
    length = 6
    def test(self, C3, C2, C1, Ccap, Cd, Cdd):
        C3_dist_O = C3.copy(atomname="O")
        C3_dist_CO = C3.copy(atomname="C")
        C1_dist = C1.copy(atomname="N")
        Ccap_dist = Ccap.copy(atomname="N")
        Cd_dist = Cd.copy(atomname="N")
        distance1a = numpdb.numdist( C3_dist_O, C1_dist )
        distance1b = numpdb.numdist( C3_dist_O, Ccap_dist )
        distance2a = numpdb.numdist( C3_dist_CO, Ccap_dist )
        distance2b = numpdb.numdist( C3_dist_CO, Cd_dist )
        return (distance1a < 3.9 or distance1b < 3.9) and (distance2a < 3.9 or distance2b < 3.9) and (C3["resname"][0]=='ASP' or C3["resname"][0]=='ASX')
    def get_dihedral(self, name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    def make_list(self, *numatoms):
        Ccap = numatoms[3]
        return (
            Ccap['resno'][0],Ccap['chain'][0],(
                self.get_dihedral( "phi", *numatoms ),
                self.get_dihedral( "psi", *numatoms )
            )
        )
class AsxTurn_type1_Test():#????????
    """
    * Residue(i)=C3 is Aspartate(D) or Asparagine(Asx)
    * H-bond between (side-chain) O of residue(i)=C3 and (main-chain) NH of residue(i+2)=C1
    * phi angles of residues C2=[-30,160]
    * psi angles of residues C3=[-140,-20], C2=[-90,40]
    Reference:
        Adel Golovin, Kim Henrick;
            MSDmotif: exploring protein sites and motifs (2008), BMC Bioinformatics 9:312
            Appendix A. Small 3D structural motifs
    """
    length = 6
    def test(self, C3, C2, C1, Ccap, Cd, Cdd):
        C3_dist_O = C3.copy(atomname="O")
        C3_dist_CO = C3.copy(atomname="C")
        C1_dist = C1.copy(atomname="N")
        Ccap_dist = Ccap.copy(atomname="N")
        Cd_dist = Cd.copy(atomname="N")
        distance1a = numpdb.numdist( C3_dist_O, C1_dist )
        distance1b = numpdb.numdist( C3_dist_O, Ccap_dist )
        distance2a = numpdb.numdist( C3_dist_CO, Ccap_dist )
        distance2b = numpdb.numdist( C3_dist_CO, Cd_dist )
        return (distance1a < 3.9 or distance1b < 3.9) and (distance2a < 3.9 or distance2b < 3.9) and (C3["resname"][0]=='SER' or C3["resname"][0]=='THR')#and C2['phi'][0]<0 and C1['phi'][0]<0 and Ccap['phi'][0]<0
    def get_dihedral(self, name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    def make_list(self, *numatoms):
        Ccap = numatoms[3]
        return (
            Ccap['resno'][0],Ccap['chain'][0],(
                self.get_dihedral( "phi", *numatoms ),
                self.get_dihedral( "psi", *numatoms )
            )
        )
class AsxTurn_type2_Test():
    """
    * Residue(i)=C3 is Aspartate(D) or Asparagine(Asx)
    * H-bond between (side-chain) O of residue(i)=C3 and (main-chain) NH of residue(i+2)=C1
    * phi angles of residues Ccap=[-99,-81], Cd=[79,91], Cdd=[-113,-87], Cddd=[-112,-88]
    * psi angles of residues C3=[-99,-81], Cd=[79,91], Cdd=[-113,-87], Cddd=[-112,-88]
    Reference:
        Adel Golovin, Kim Henrick;
            MSDmotif: exploring protein sites and motifs (2008), BMC Bioinformatics 9:312
            Appendix A. Small 3D structural motifs
    """
    length = 6
    def test(self, C3, C2, C1, Ccap, Cd, Cdd):
        C3_dist_O = C3.copy(atomname="O")
        C3_dist_CO = C3.copy(atomname="C")
        C1_dist = C1.copy(atomname="N")
        Ccap_dist = Ccap.copy(atomname="N")
        Cd_dist = Cd.copy(atomname="N")
        distance1a = numpdb.numdist( C3_dist_O, C1_dist )
        distance1b = numpdb.numdist( C3_dist_O, Ccap_dist )
        distance2a = numpdb.numdist( C3_dist_CO, Ccap_dist )
        distance2b = numpdb.numdist( C3_dist_CO, Cd_dist )
        return (distance1a < 3.9 or distance1b < 3.9) and (distance2a < 3.9 or distance2b < 3.9) and (C3["resname"][0]=='SER' or C3["resname"][0]=='THR')#and C2['phi'][0]<0 and C1['phi'][0]<0 and Ccap['phi'][0]<0
    def get_dihedral(self, name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    def make_list(self, *numatoms):
        Ccap = numatoms[3]
        return (
            Ccap['resno'][0],Ccap['chain'][0],(
                self.get_dihedral( "phi", *numatoms ),
                self.get_dihedral( "psi", *numatoms )
            )
        )
class STTest():
    """
    * Residue(i)=C3 is Serine(S) or Threonine(T)
    * H-bond between (side-chain) O of residue(i)=C3 and (main-chain) NH of residue(i+2)=C1 or residue(i+3)=Ccap
    * H-bond between CO of residue(i)=C3 and NH of residue(i+3)=Ccap or residue(i+4)=Cd
    Reference:
        Adel Golovin, Kim Henrick;
            MSDmotif: exploring protein sites and motifs (2008), BMC Bioinformatics 9:312
            Appendix A. Small 3D structural motifs
    """
    length = 6
    def test(self, C3, C2, C1, Ccap, Cd, Cdd):
        C3_dist_O = C3.copy(atomname="O")
        C3_dist_CO = C3.copy(atomname="C")
        C1_dist = C1.copy(atomname="N")
        Ccap_dist = Ccap.copy(atomname="N")
        Cd_dist = Cd.copy(atomname="N")
        distance1a = numpdb.numdist( C3_dist_O, C1_dist )
        distance1b = numpdb.numdist( C3_dist_O, Ccap_dist )
        distance2a = numpdb.numdist( C3_dist_CO, Ccap_dist )
        distance2b = numpdb.numdist( C3_dist_CO, Cd_dist )
        return (distance1a < 3.9 or distance1b < 3.9) and (distance2a < 3.9 or distance2b < 3.9) and (C3["resname"][0]=='SER' or C3["resname"][0]=='THR')#and C2['phi'][0]<0 and C1['phi'][0]<0 and Ccap['phi'][0]<0
    def get_dihedral(self, name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    def make_list(self, *numatoms):
        Ccap = numatoms[3]
        return (
            Ccap['resno'][0],Ccap['chain'][0],(
                self.get_dihedral( "phi", *numatoms ),
                self.get_dihedral( "psi", *numatoms )
            )
        )
    
    



class AlphaLTest():
    """
    * H-bond between CO of residue C3 and NH of residue Cd
    * phi angles of residues Ccap=[-99,-81], Cd=[79,91], Cdd=[-113,-87], Cddd=[-112,-88]
    * psi angles of residues Ccap=[-16,-2], Cd=[3,21], Cdd=[110,144]
    Reference:
        R Aurora, GD Rose;
            Helix capping (1998), Protein Sci. 7(1):21-38
    """
    length = 6
    def test(self, C3, C2, C1, Ccap, Cd, Cdd):#, Cddd):
        C3_dist = C3.copy(atomname="O")
        Cd_dist = Cd.copy(atomname="N")
        distance1 = numpdb.numdist( C3_dist, Cd_dist )
        return distance1<3.9 and -99<Ccap['phi'][0]<-81 and 79<Cd['phi'][0]<91 and -113<Cdd['phi'][0]<-87 and -16<Ccap['psi'][0]<-2 and 3<Cd['psi'][0]<21 and 110<Cdd['psi'][0]<144# and -112<Cddd['phi'][0]<-88 
    def get_dihedral(self, name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    def make_list(self, *numatoms):
        Ccap = numatoms[3]
        return (
            Ccap['resno'][0],Ccap['chain'][0],(
                self.get_dihedral( "phi", *numatoms ),
                self.get_dihedral( "psi", *numatoms )
            )
        )

class AlphaLTest_fitted():
    """
    * H-bond between CO of residue C3 and NH of residue Cd
    * phi angles of residues C3=[-100,50], C2=[-100,50], C1=[-100,50], Ccap=[-150,-25], Cd=[0,180]
    * psi angles of residues C3=[-100,50], C2=[-100,50], C1=[-100,50], Ccap=[-75,50]
    Reference:
        self
    """
    length = 8
    def test(self, C5, C4, C3, C2, C1, Ccap, Cd, Cdd):#, Cddd):
        C3_dist = C3.copy(atomname="O")
        C2_dist = C2.copy(atomname="O")
        Cd_dist = Cd.copy(atomname="N")
        distance1 = numpdb.numdist( C3_dist, Cd_dist )
        distance2 = numpdb.numdist( C2_dist, Cd_dist )
        return distance1<3.4 and distance2<3.4 and -100<C3['phi'][0]<50 and C4["sstruc"][0]=="H" and C5["sstruc"][0]=="H" \
            and -100<C3['psi'][0]<50  and -100<C2['phi'][0]<50 and C3["sstruc"][0]=="H" \
            and -100<C2['psi'][0]<50  and -100<C1['phi'][0]<50 and C2["sstruc"][0]=="H"  \
            and -100<C1['psi'][0]<50  and -150<Ccap['phi'][0]<-25 \
            and -75<Ccap['psi'][0]<50  and 0<Cd['phi'][0]<180  #and (-180<Cd['psi'][0]<-100 or 100<Cd['psi'][0]<180)
    def get_dihedral(self, name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    def make_list(self, *numatoms):
        Ccap = numatoms[5]
        return (
            Ccap['resno'][0],Ccap['chain'][0],(
                self.get_dihedral( "phi", *numatoms )[2:],
                self.get_dihedral( "psi", *numatoms )[2:]
            )
        )






npdb_dict = {}
def get_npdb( pdb_id ):
    pdb_repos='/home/student/Johanna/Projekte/caps/cap-referenceset/'
    if pdb_id not in npdb_dict:
        path = pdb_repos+pdb_id+'.pdb'
        npdb_dict[pdb_id] = numpdb.NumPdb( path, features={
            "phi_psi": True,
            "sstruc": True
        })
    return npdb_dict[pdb_id]



def get_refset( pdbid_list, outdir ):
    pdbid=[]
    with open(pdbid_list) as pdbids:
        for line in pdbids:
            if not line.startswith('#'):
                pdbid.append( line.upper()[0:4] )
    refdir = os.path.join(outdir, 'refdir')
    try: os.makedirs(refdir)
    except: pass
    for i in pdbid:
        pdb_download(i, os.path.join(refdir, i+'.pdb'))
    return pdbid
from opm import Opm
import opm
from utils.math import point_plane_dist


def get_TM_regions( pdbid, min_helix_length, outpath ):
    opm = Opm( pdbid, output_dir=outpath+"/opm/" )
    mplanes = opm.get_planes()
    npdb = NumPdb( opm.processed_file )
    dist = abs( point_plane_dist( mplanes[1][0], mplanes[0] ) )
    npdb_tm_dict = {}
    npdb_sol_dict = {}
    sele = npdb.sele()
    i = 0; j = 0; k = 0
    tm_end = []; tm_start = []
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
        if flag:#solvent
            if j>=min_helix_length:
                tm_end.append((numa["chain"][0],numa["resno"][0]))
            tempstr = ((numa["chain"][0],numa["resno"][0]))
            j = 0; k = 0
        else:#tm
            j += 1
            k += 1
            if k>=min_helix_length:
                tm_start.append(tempstr)
                k = -99
    tm_regions = zip(tm_start, tm_end)
    tm_regions.write( outpath+ "tm_regions.txt" )
    npdb_sol_dict[ inpdb[0][0] ] = npdb.copy( sele=sele )
    npdb.copy( sele=sele ).write( outpath+ "sol_region.pdb" )
    np.logical_not( sele, sele )
    npdb_tm_dict[ inpdb[0][0] ] = npdb.copy( sele=sele )
    npdb.copy( sele=sele ).write( outpath+ "tm_region.pdb"  )
    return npdb_tm_dict, npdb_sol_dict, tm_regions

def get_helix_regions( pdbid, refdir, min_helix_length=6, dssp=False ):
    helix_regions = []
    def get_file_path( pdbid ):
        return refdir+pdbid+'.pdb'
    if dssp:
        numa = get_npdb( pdbid );
        curresno=0; count=0; startresno=0
        begin=[];end=[]
        for atom in numa:
            if atom['sstruc']=='H':
                if curresno!=atom["resno"]:
                    curresno = atom["resno"]
                    if count==0: startresno=atom["resno"]
                    count += 1
            else:
                if count >= min_helix_length:
                    begin.append( ( atom["chain"],startresno ) )
                    end.append( ( atom["chain"],curresno ) )
                count=0
        pdbidd = []; i = 0
        while i<len(begin):
            pdbidd.append(pdbid)
            i +=1
        helix_regions = zip(pdbidd,zip(begin, end))#under construction
    else:#helixenden aus der pdb
        with open(get_file_path(pdbid), 'r') as fp:
            for line in fp:
                if line.startswith('HELIX'):
                    # chain, resno, resname
                    start = ( line[19], try_int( line[21:25] ), line[15:18] )
                    end = ( line[31], try_int( line[33:37] ), line[27:30] )
                    if start>end: start, end = end, start
                    helix_regions.append( ( pdbid,start,end ) )
                elif line.startswith('SHEET'):
                    break;
    return helix_regions
#print get_helix_regions('1ABA','/home/student/Johanna/Projekte/caps/cap-referenceset/')

def calc_distances( pdbid, chain, refresno, dist_range, side_chain=False ):
    def distnum( resno, atomname ):
        return numa.copy( chain=chain, resno=resno, atomname=atomname )
    dist_list = []
    numa = get_npdb( pdbid )
    if side_chain:
        pass
    else:
        for resno in range( refresno-dist_range, refresno+dist_range ):
            if resno <= 0:
                dist_list.append( ( pdbid+chain+str(resno), None ) )
            else:
                distance = numpdb.numdist( distnum( resno, "O" ), distnum( refresno, "N" ) )
                distance2 = numpdb.numdist( distnum( resno, "N" ), distnum( refresno, "O" ) )
                if distance <= 3.5:
                    dist_list.append( ( pdbid+chain+str(resno), distance ) )
                elif distance2 <= 3.5:
                    dist_list.append( ( pdbid+chain+str(resno), distance2 ) )
                else: dist_list.append( ( pdbid+chain+str(resno), None ) )
    return dist_list

def superpose2( refpdbid, refchain, refresno, pdbid, chain, resno, superpose_range=1 ):#superpose_range = out of helix
    numa = get_npdb( refpdbid ); numa2 = get_npdb( pdbid )
    sele_ref = {
            "chain": refchain, "resno": [ refresno-4, refresno+superpose_range ]
        }
    sele = {
            "chain": chain, "resno": [ resno-4, resno+superpose_range ]
        }
    #print sele_ref, sele
    try:
        rmsd = numpdb.superpose(
            numa, numa2, sele, sele_ref, align=False#, verbose=False
        )
    except Exception:
        rmsd=None
    
    try:
        return rmsd.rmsd
    except Exception: #to do: length differ, cannot superpose
        return 999
    
#print superpose2( '1ABA', 'A', 15, '1ABA', 'A', 30 )
def get_angles( pdbid, chain, resno ):
    def distnum( chain,resno ):
        return numa.copy( chain=chain, resno=resno )
    def get_dihedral(name, *numatoms):
        return np.hstack( map( lambda x: x[name][0], numatoms ) )
    numa = get_npdb( pdbid );
    phi = get_dihedral( "phi", distnum(chain, resno) )
    psi = get_dihedral( "psi", distnum(chain, resno) )
    return ( phi, psi )
#print get_angles('1ABA', 'A', 15)

def get_infos( ccap_range, ccap, dist_range ):
    pdbid, chain, resno = ccap
    dist_list = []; angle_list = []
    for i in range(resno-ccap_range, resno+ccap_range+1):
        if i <= 0:
            dist_list.append( ( i,None ) )
            angle_lsit.append( ( None,None ) )
        else:
            dist_list.append( calc_distances( pdbid, chain, i,  dist_range ) )
            angle_list.append( get_angles( pdbid, chain, i ) )
    #todo: save lists to file
    return dist_list, angle_list



outdir = '/home/student/Johanna/Projekte/caps/tm-cluster/'
pdbid_list = os.path.join(outdir,'pdb-list.txt')
#print get_refset(pdbid_list, outdir)
#print get_TM_regions( '2GIF', 6, '/home/student/Johanna/Projekte/caps/tm-cluster/')






#clusterbysuperpose
def cluster_by_superpose():
    #inputs
    outputdir='/home/student/Johanna/Projekte/caps/cluster_by_superpose/'
    refdir='/home/student/Johanna/Projekte/caps/cap-referenceset/'
    min_helix_length=8
    
    
    all_ccaps = []
    all_pdbs = [i[0:-4] for i in os.listdir(refdir)]
    for pdbid in all_pdbs:
        #get ccap:
        helix_regions = get_helix_regions( pdbid, refdir, min_helix_length)
        for elem in helix_regions:
            all_ccaps.append( ( elem[0],elem[1][0],elem[1][1] ) )
    all_ccaps.sort()
    superpose_matrix = np.empty([len(all_ccaps),len(all_ccaps)])
    #superpose_list = []
    for index, refpdbids in enumerate(all_ccaps):
        refpdbid, refchain, refresno = refpdbids
        for index2, pdbids in enumerate(all_ccaps):
            if index2>index:
                pdbid, chain, resno = pdbids
                print index, len(all_ccaps), index2
                rmsd = superpose2(refpdbid, refchain, refresno, pdbid, chain, resno, superpose_range=1)
                #superpose_list.append( ( rmsd, ( re ), (  ) ) )
                superpose_matrix[index][index2]=rmsd
                superpose_matrix[index2][index]=rmsd
            elif index2==index:
                superpose_matrix[index][index2]=999
    curr_list=[]
    # todo: cluster!!! spalte/zeile mit niedrigstem Wert im bereich [0:0.2?]
    # wenn gleich: gesamte spalte/zeile
    # die raus, die in dem Bereich drin sind
    # neue matrix schreiben, an all_ccaps denken!
    print np.amin(superpose_matrix)
    print np.argmin(superpose_matrix)
    print superpose_matrix
    
    
    pass
cluster_by_superpose()
def save_newcluster( newcluster, refmotif, refdir, outputdir, overlap ):
    def _make_jspt_color_motif( pdbid, chain, resno, jspt_file ):
        with open( jspt_file, "w" ) as fp:
            s = "select {%s:%s}; color @{color('roygb', 0, %s, %s)}; wireframe 0.3; \n" % (resno, chain, 0,0)
            fp.write( s )
            fp.write( "select none; background black;" )
    def _make_provi_file( pbdid, refdir, provi_file, jspt_file ):
        with open( provi_file, "w" ) as fp:
            s = '[ {"filename": "%s"},{"filename":  "%s"}]' % ('../../../cap-referenceset/'+pdbid+'.pdb', jspt_file)
            fp.write( s )
    
    if overlap: outdir2=outputdir+'cluster-overlap/'+refmotif+'/'
    else: outdir2=outputdir+'cluster-no-overlap/'+refmotif+'/'
    try: os.makedirs(outdir2)
    except: pass
    for motif in newcluster:
        pdbid=motif[0:4]; chain=motif[4]; resno=motif[5:]
        jspt_file=outdir2+pdbid+chain+resno+'_color_motif.jspt'
        provi_file=outdir2+pdbid+chain+resno+'_motif.provi'
        _make_jspt_color_motif(pdbid, chain, resno, jspt_file)
        _make_provi_file( pdbid, refdir, provi_file, pdbid+chain+resno+'_color_motif.jspt' )
def make_clusters_with_overlap(  ):
    #all_motifs=get_all_motifs( indir )
    motifs_file='/home/student/Johanna/Projekte/caps/norange/2013-08-07_14-59-04_all_motifs.txt'
    outputdir='/home/student/Johanna/Projekte/caps/cluster/'
    refdir='/home/student/Johanna/Projekte/caps/cap-referenceset/'
    all_motifs = []
    with open(motifs_file, 'r') as fo:
        for motifline in fo:
            pdb_id1, chain1, cd1 = motifline.split()
            all_motifs.append(pdb_id1+chain1+str(cd1))
    curr_motif_list=[]
    curr_motif_list=all_motifs
    
    while curr_motif_list!=[]:
        curr_treshold=get_treshold(curr_motif_list, '0,0.2')
        max_treshold=[]
        j=0
        while j<=len(curr_treshold)-1 and (not curr_treshold[j][0]<curr_treshold[0][0]):
            max_treshold.append(curr_treshold[j])
            j += 1
        newcluster_no_overlap = superpose(sorted(max_treshold, key=lambda tup: tup[2])[0][1], curr_motif_list, '0,0.2' )
        save_newcluster( newcluster_no_overlap, sorted(max_treshold, key=lambda tup: tup[2])[0][1], refdir, outputdir, False )
        for i in max_treshold:
            newcluster_overlap = superpose(i[1], all_motifs, '0,0.2' )
            print i
            save_newcluster( newcluster_overlap, i[1], refdir, outputdir, True )
            curr_motif_list=list((set(curr_motif_list)).difference((set(curr_motif_list)).intersection(set(newcluster_overlap))))


jo=kojo
print jo
#clusterbydistances - cluster by superpose























motifs_dict = {
    "alpha_beta_motif": AlphaBetaMotifTest(),
    "alpha_L_motif": AlphaLTest(),
    "alpha_L_motif-fitted": AlphaLTest_fitted(),
    "st_motif": STTest()
}
npdb_dict = {}
def get_npdb( pdb_id ):
    pdb_repos='/home/student/Johanna/Projekte/caps/cap-referenceset/'
    if pdb_id not in npdb_dict:
        path = pdb_repos+pdb_id+'.pdb'
        npdb_dict[pdb_id] = numpdb.NumPdb( path, features={
            "phi_psi": False,
            "sstruc": False
        })
    return npdb_dict[pdb_id]
def _superpose_all( numa, sele_motif, numa2, sele_motif2):
    
    sele_ref = {
            "chain": sele_motif2['chain'], "resno": [ sele_motif2['cd']-4, sele_motif2['cd'] ]
        }
    sele = {
            "chain": sele_motif['chain'], "resno": [ sele_motif['cd']-4, sele_motif['cd'] ]
        }
    old_stdout = sys.stdout
    sys.stdout = stdout = StringIO()
    numpdb.superpose(
            numa, numa2, sele, sele_ref, align=False
        )
    sys.stdout = old_stdout
    rmsd = stdout.getvalue().split()[1]
    sys.stdout = sys.stdout
    return float(rmsd)
def get_treshold( curr_motif_list, t='0,0.2' ):
    motifs = curr_motif_list
    arr= np.empty( (len(motifs)*(len(motifs)-1))/2 )
    k=0; j=0
    curr_treshold=[]
    for i, mo in enumerate(motifs):
        pdb_id1=mo[0:4]; chain1=mo[4]; cd1=mo[5:]
        n1=get_npdb(pdb_id1)
        al=[]; kj=0
        while kj<len(motifs):
            pdb_id2=motifs[kj][0:4]; chain2=motifs[kj][4]; cd2=motifs[kj][5:]
            n2=get_npdb(pdb_id2)
            rmsd=_superpose_all(
                    n1, {'chain':chain1, 'cd':int(cd1)},
                    n2 , {'chain':chain2, 'cd':int(cd2) }
                )
            kj +=1
            al.append(rmsd)
        hist2db, phiedgesb, psiedgesb = np.histogram2d( al, al, bins=1, range=[[t.split(',')[0],t.split(',')[1]],[t.split(',')[0],t.split(',')[1]]])
        list_histb=[]
        for i, line in enumerate(hist2db):
            list_histb.append(hist2db[i][i])
        curr_treshold.append([list_histb[0],mo, sum(al)/len(al)])
        
    return sorted(curr_treshold, reverse=True)

def superpose( ref, all_motifs, t='0,0.2' ):
    
    numa_ref=get_npdb( ref[0:4] )
    chain_ref=ref[4]
    cd_ref=ref[5:]
    sele_ref = {
            "chain": chain_ref, "resno": [ int(cd_ref)-4, int(cd_ref) ]
        }
    newcluster=[]
    for motif in all_motifs:
        numa=get_npdb( motif[0:4] )
        chain=motif[4]
        cd=motif[5:]
        sele = {
                "chain": chain, "resno": [ int(cd)-4, int(cd) ]
            }
        old_stdout = sys.stdout
        sys.stdout = stdout = StringIO()
        numpdb.superpose(
                numa, numa_ref, sele, sele_ref, align=False
            )
        sys.stdout = old_stdout
        rmsd = stdout.getvalue().split()[1]
        sys.stdout = sys.stdout
        if float(rmsd) < float(t.split(',')[1]):
            newcluster.append(motif)
    return sorted(newcluster)

def save_newcluster( newcluster, refmotif, refdir, outputdir, overlap ):
    def _make_jspt_color_motif( pdbid, chain, resno, jspt_file ):
        with open( jspt_file, "w" ) as fp:
            s = "select {%s:%s}; color @{color('roygb', 0, %s, %s)}; wireframe 0.3; \n" % (resno, chain, 0,0)
            fp.write( s )
            fp.write( "select none; background black;" )
    def _make_provi_file( pbdid, refdir, provi_file, jspt_file ):
        with open( provi_file, "w" ) as fp:
            s = '[ {"filename": "%s"},{"filename":  "%s"}]' % ('../../../cap-referenceset/'+pdbid+'.pdb', jspt_file)
            fp.write( s )
    
    if overlap: outdir2=outputdir+'cluster-overlap/'+refmotif+'/'
    else: outdir2=outputdir+'cluster-no-overlap/'+refmotif+'/'
    try: os.makedirs(outdir2)
    except: pass
    for motif in newcluster:
        pdbid=motif[0:4]; chain=motif[4]; resno=motif[5:]
        jspt_file=outdir2+pdbid+chain+resno+'_color_motif.jspt'
        provi_file=outdir2+pdbid+chain+resno+'_motif.provi'
        _make_jspt_color_motif(pdbid, chain, resno, jspt_file)
        _make_provi_file( pdbid, refdir, provi_file, pdbid+chain+resno+'_color_motif.jspt' )
        
    



def make_clusters_with_overlap(  ):
    #all_motifs=get_all_motifs( indir )
    motifs_file='/home/student/Johanna/Projekte/caps/norange/2013-08-07_14-59-04_all_motifs.txt'
    outputdir='/home/student/Johanna/Projekte/caps/cluster/'
    refdir='/home/student/Johanna/Projekte/caps/cap-referenceset/'
    all_motifs = []
    with open(motifs_file, 'r') as fo:
        for motifline in fo:
            pdb_id1, chain1, cd1 = motifline.split()
            all_motifs.append(pdb_id1+chain1+str(cd1))
    curr_motif_list=[]
    curr_motif_list=all_motifs
    
    while curr_motif_list!=[]:
        curr_treshold=get_treshold(curr_motif_list, '0,0.2')
        max_treshold=[]
        j=0
        while j<=len(curr_treshold)-1 and (not curr_treshold[j][0]<curr_treshold[0][0]):
            max_treshold.append(curr_treshold[j])
            j += 1
        newcluster_no_overlap = superpose(sorted(max_treshold, key=lambda tup: tup[2])[0][1], curr_motif_list, '0,0.2' )
        save_newcluster( newcluster_no_overlap, sorted(max_treshold, key=lambda tup: tup[2])[0][1], refdir, outputdir, False )
        for i in max_treshold:
            newcluster_overlap = superpose(i[1], all_motifs, '0,0.2' )
            print i
            save_newcluster( newcluster_overlap, i[1], refdir, outputdir, True )
            curr_motif_list=list((set(curr_motif_list)).difference((set(curr_motif_list)).intersection(set(newcluster_overlap))))





def seco_cluster():
    coords = np.array([ [ d['x'], d['y'], d['z'] ] for d in vert ])
    fd = scipy.cluster.hierarchy.fclusterdata(
        coords, pr, criterion='distance',
        metric='euclidean', method='average'
    )
    clust = collections.defaultdict( list )
    for i, x in enumerate(fd):
        clust[x].append( coords[i] )    
    pass




def _superpose_test( numa, sele_motif ):

    numa_ref=numpdb.NumPdb("/home/student/Johanna/Projekte/caps/cap-referenceset/3COX.pdb")
    chain=sele_motif['chain']
    cd=sele_motif['cd']
    sele_ref = {
            "chain": "A", "resno": [ 304-4, 304 ]
        }
    sele = {
            "chain": chain, "resno": [ cd-4, cd ]
        }
    old_stdout = sys.stdout
    sys.stdout = stdout = StringIO()
    numpdb.superpose(
            numa, numa_ref, sele, sele_ref, align=False
        )
    sys.stdout = old_stdout
    rmsd = stdout.getvalue().split()[1]
    sys.stdout = sys.stdout
    return float(rmsd)

def _superpose_all( numa, sele_motif, numa2, sele_motif2):
    
    sele_ref = {
            "chain": sele_motif2['chain'], "resno": [ sele_motif2['cd']-4, sele_motif2['cd'] ]
        }
    sele = {
            "chain": sele_motif['chain'], "resno": [ sele_motif['cd']-4, sele_motif['cd'] ]
        }
    old_stdout = sys.stdout
    sys.stdout = stdout = StringIO()
    numpdb.superpose(
            numa, numa2, sele, sele_ref, align=False
        )
    sys.stdout = old_stdout
    rmsd = stdout.getvalue().split()[1]
    sys.stdout = sys.stdout
    return float(rmsd)



def supallall( motifs_file, pdb_repos, rmsdrange, *args ):
    npdb_dict = {}
    def get_npdb( pdb_id ):
        if pdb_id not in npdb_dict:
            path = pdb_repos+pdb_id+'.pdb'
            npdb_dict[pdb_id] = numpdb.NumPdb( path, features={
                "phi_psi": False,
                "sstruc": False
            })
        return npdb_dict[pdb_id]
    motifs = []
    with open(motifs_file, 'r') as fo:
        for motifline in fo:
            pdb_id1, chain1, cd1 = motifline.split()
            motifs.append((pdb_id1, chain1, cd1))    
    arr= np.empty( (len(motifs)*(len(motifs)-1))/2 )
    k=0; highrmsd=0; j=0
    array_all={}
    for i, mo in enumerate(motifs):
        pdb_id1, chain1, cd1 = mo
        n1=get_npdb(pdb_id1)
        j=i+1
        
        while j<len(motifs):
            pdb_id2, chain2, cd2 = motifs[j]
            n2=get_npdb(pdb_id2)
            rmsd=_superpose_all(
                    n1, {'chain':chain1, 'cd':int(cd1)},
                    n2 , {'chain':chain2, 'cd':int(cd2) }
                )

            if rmsd>highrmsd:
                highrmsd=rmsd
            arr[k]=rmsd
            j += 1
            k +=1
        al=[]; kj=0
        while kj<len(motifs):
            pdb_id2, chain2, cd2 = motifs[kj]
            n2=get_npdb(pdb_id2)
            rmsd=_superpose_all(
                    n1, {'chain':chain1, 'cd':int(cd1)},
                    n2 , {'chain':chain2, 'cd':int(cd2) }
                )
            kj +=1
            al.append(rmsd)
        hist2db, phiedgesb, psiedgesb = np.histogram2d( al, al, bins=1, range=[[0,0.2],[0,0.2]])
        list_histb=[]
        for i, line in enumerate(hist2db):
            list_histb.append(hist2db[i][i])
        array_all[str(pdb_id1 + chain1 + cd1)]=list_histb
    for dir in args:a=dir
    name2=a+"superpose"+'_'+strftime("%Y-%m-%d_%H-%M-%S", gmtime())+'.txt'
    with open( name2 , "w" ) as fp2:
        for j in phiedgesb:
            fp2.write(str(j)+'\t')
        fp2.write('\n')
        for el in array_all:
            for i in array_all[el]:
                fp2.write(str(int(i))+'\t')
            
            fp2.write(el +'\n')
        
    if rmsdrange[1]==0: rmsdrange[0]=0; rmsdrange[1]=highrmsd
    hist2d, phiedges, psiedges = np.histogram2d( arr, arr, bins=25,
                                                range=[[rmsdrange[0],rmsdrange[1]],[rmsdrange[0],rmsdrange[1]]]
                                                #[[0,highrmsd+0.1],[0,highrmsd+0.1]])
                                                )
    list_hist=[]
    for i, line in enumerate(hist2d):
        list_hist.append(hist2d[i][i])
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ind = np.arange(len(list_hist))
    ax.bar(ind, list_hist)
    ax.set_xticks(ind+0.75)
    ax.set_xticklabels( phiedges, fontsize=8 )
    plt.ylabel('#motifs')
    plt.xlabel('RMSD')
    fig.autofmt_xdate()
    a=""
    for dir in args:a=dir
    name=a+"superpose"+'_'+strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    plt.savefig(name, dpi=300)
    pass



def _make_record( motif, name, motif_line, i, rmsd, l2f_str, m2f_str ):
    return MotifRecord(
        motif, name, motif_line[1], motif_line[0], rmsd,
        list(motif_line[2][0]), list(motif_line[2][1]), l2f_str, m2f_str, i
    )

def _find_all_motifs( npdb ):
    # takes a numpy-pdb to the motif-classes and
    # gets a dicc with the motif-name and the output of the motif-classes
    motif_data = defaultdict(list); records = []
    for name, motif_obj in motifs_dict.iteritems():
        motif_data[ name ] = []
        for numatoms in npdb.iter_resno2( motif_obj.length ):
            if motif_obj.test(*numatoms):
                motif_data[ name ].append( motif_obj.make_list(*numatoms) )
    return motif_data

def find_motifs( infile, motif_type, m2f, l2f, rmsdrange):
    # generate a numpy-pdb
    # get the motif-dicc and parses the motifs to Motif-Record
    records = []
    npdb=numpdb.NumPdb(infile)
    pdb_id = utils.path.stem( infile )
    for motif in motif_type:
        motif_found=_find_all_motifs(npdb)[motif]
        if motif_found!=[]:
            for i, elem in enumerate(motif_found):
                rmsd=_superpose_test( npdb, {'chain':elem[1], 'cd':elem[0]} )
                resstr=""
                l2f_str=''
                if l2f=='True':
                    for i in range(elem[0]-4, elem[0]):
                        resstr=resstr+bio.AA1.get(npdb.iter_resno2(2)[i-4][0][0]["resname"], "X")
                    l2f_str=">"+pdb_id+"|"+elem[1]+str(elem[0])+"\n"+resstr+'\n'
                m2f_str=''
                if m2f=='True':
                    m2f_str=pdb_id+" "+elem[1]+" "+str(elem[0])+"\n"
                if rmsdrange[0]!='':
                    if rmsd > float(rmsdrange[0]) and rmsd < float(rmsdrange[1]):
                        records.append(_make_record(motif, pdb_id, elem, i, rmsd, l2f_str, m2f_str))
                else: records.append(_make_record(motif, pdb_id, elem, i, rmsd, l2f_str, m2f_str))
    return records
#Muscle('/home/student/Johanna/Projekte/caps/fasta3.fasta')

def _motif_out_record( records ):
    motifdic=defaultdict(list)
    for r in records:
        motifdic[r.motif].append(r)
    return motifdic

def plot_found_motifs(records, plot, motif, outdir):
    ra=_motif_out_record(records)
    for mo in motif:
        _overplot( ra[mo], plot, mo, outdir)

def _overplot( records, plot, motife, outdir ):
    n = 0; m = 0; axes=None
    if plot=="big": fig, axes  = plt.subplots( nrows=2, ncols=3, figsize=(15,15) )
    for index, pos in enumerate(['C3','C2','C1','CCap','C\'','C\'\'']):
        phi=[]; psi=[];
        for motif_pos in records:
            phi=np.hstack((phi,motif_pos.phi_angles[index]))
            psi=np.hstack((psi,motif_pos.psi_angles[index]))
        if not (phi==[] or psi==[]):
            sele=(np.isnan(phi)==False) & (np.isnan(psi)==False)
            phi = phi[sele]
            psi = psi[sele]
            _rama_plot( phi, psi, pos,motife,plot,n, m,axes, outdir )
            m=m+1
            if m==3:
                m=0
                n=n+1
    if plot=="big":
        fig.tight_layout()
        name=outdir+motife+'-motif_all'+'_'+strftime("%Y-%m-%d_%H-%M-%S", gmtime())
        plt.suptitle('Ramachandrian plot for '+motife+' motif')
        plt.savefig(name, dpi=300)
        plt.close()
    pass

def _rama_plot( phi, psi, pos, motif, plot, n, m, axes, *args ):
    hist2d, phiedges, psiedges = np.histogram2d( phi, psi, bins=180, range=[[-180,180],[-180,180]] )
    extent=[-180,180,-180,180]
    if plot=="big":
        axes[n,m].imshow( np.flipud(hist2d), extent=extent, interpolation='nearest', cmap=cm.Reds)
        axes[n,m].set_title('residue '+pos)
        axes[n,m].set_ylabel('phi angles')
        axes[n,m].set_xlabel('psi angles')
    else:
        plt.imshow( np.flipud(hist2d), extent=extent, interpolation='nearest', cmap=cm.Reds)
        plt.colorbar()
        plt.title('Ramachandrian plot for '+motif+' motif at residue '+pos)
        plt.ylabel('phi angles')
        plt.xlabel('psi angles')
        plt.plotting()
        for dir in args:a=dir
        name=a+motif+pos+'_'+strftime("%Y-%m-%d_%H-%M-%S", gmtime())
        plt.savefig(name, dpi=300)
        plt.close()

def _make_struc_record( name, info_line, i ):
    #    'pdb_id', 'chain', 'resno', 'phi_angles', 'psi_angles', 'distances', 'no'
    return StrucRecord(
        name, info_line[0], info_line[1], 
        list(info_line[2]), list(info_line[3]),
        list(info_line[4]), i
    )

def _motif_out_struc_record( records ):
    infodic=defaultdict(list)
    for r in records:
        infodic[r.resno].append(r)
    return infodic

def get_info( infile, residue ):
    # generate a numpy-pdb
    # get the motif-dicc and parses the motifs to Motif-Record
    infos = []
    npdb=numpdb.NumPdb(infile)
    pdb_id = utils.path.stem( infile )

    for i, res in enumerate(residue):
        distances=[]; phi_angles=[]; psi_angles=[]
        chain = res.split(':')[0]
        resi = int(res.split(':')[1])
        
        for numatoms in npdb.iter_resno2( 7 ):
            if numatoms[3]['resno'][0]==resi and numatoms[3]['chain'][0]==chain:
                for el in numatoms:
                    phi_angles.append( el[0]['phi'] )
                    psi_angles.append( el[0]['psi'] )
                    dist = numpdb.numdist( el.copy(atomname='O'), numatoms[4].copy(atomname='N') )
                    if 0 < dist < 13.9:
                        distances.append((el[0]['resno'], numatoms[4]['resno'][0], dist))
                break;
        elem=[chain, resi, phi_angles, psi_angles, distances]
        infos.append(_make_struc_record(pdb_id, elem, i))
    return infos

def _histo_plot( dist_list, motif, *args ):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    y = [np.mean(dist_list[0]), np.mean(dist_list[1]), np.mean(dist_list[2]), np.mean(dist_list[3])]
    N = len(y)
    ind = range(N)
    err = [np.std(dist_list[0])/np.sqrt(len(dist_list[0])), np.std(dist_list[1])/np.sqrt(len(dist_list[1])), np.std(dist_list[2])/np.sqrt(len(dist_list[2])), np.std(dist_list[3])/np.sqrt(len(dist_list[3]))]
    # See note below on the breakdown of this command
    ax.bar(ind, y, facecolor='#777777',
           align='center', yerr=err, ecolor='black')
    ax.set_ylabel('Distance')
    ax.set_title('Distances from CO to N',fontstyle='italic')
    ax.set_xticks(ind)
    group_labels = ['C3 to C\'', 'C2 to C\'',
                     'C1 to C\'', 'Ccap to C\'']
    ax.set_xticklabels(group_labels)
    fig.autofmt_xdate()
    plt.plotting()
    for dir in args:a=dir
    name=a+motif+'_Distance_'+strftime("%Y-%m-%d_%H-%M-%S", gmtime())
    plt.savefig(name, dpi=300)
    plt.close()

def ploti( outdir, motife, plot ):
    phi_list=[]; psi_list=[]; distance_list=[]
    for dir, subdir, files in os.walk(outdir):
        for fi in files:
            if fi[-4:]=='json' and not fi=='structuregetter.json':
                with open(fi) as json_data:
                    data = json.loads(json_data.read())
                    phi_list.append(data[0]["phi_angles"])
                    psi_list.append(data[0]["psi_angles"])
                    distance_list.append(data[0]["distances"])
    n = 0; m = 0; axes=None
    if plot=="big": fig, axes  = plt.subplots( nrows=3, ncols=3, figsize=(15,15) )
    for index, pos in enumerate(['C3','C2','C1','CCap','C\'','C\'\'','C\'\'\'']):
        phi=[]; psi=[];
        for motif_pos in phi_list:
            phi=np.hstack((phi,motif_pos[index]))
        for motif_pos in psi_list:
            psi=np.hstack((psi,motif_pos[index]))
        if not (phi==[] or psi==[]):
            sele=(np.isnan(phi)==False) & (np.isnan(psi)==False)
            phi = phi[sele]
            psi = psi[sele]
            _rama_plot( phi, psi, pos,motife,plot,n, m,axes, outdir )
            m=m+1
            if m==3:
                m=0
                n=n+1
    if plot=="big":
        fig.tight_layout()
        name=outdir+motife+'-motif_all'+'_'+strftime("%Y-%m-%d_%H-%M-%S", gmtime())
        plt.suptitle('Ramachandrian plot for '+motife+' motif')
        plt.savefig(name, dpi=300)
        plt.close()
    #distance-plot
    dist_list=[]
    for index, pos in enumerate(['C3','C2','C1','CCap']):
        dist=[];
        for motif_pos in distance_list:
            dist=np.hstack((dist,motif_pos[index][2]))
        dist_list.append(dist)
    if not dist_list==[]:
        _histo_plot( dist_list, motife, outdir )
    pass


class CapsMotifFinder( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """A tool that finds cap motifs. """
    args = [
        _( "inputs", type="text" ,
            help="input: file, files (seperated with a ',') or one directory, "
            "e.g. '~/Downloads' or '~/Downloads/3DQB.pdb,~/Downloads/3SN6'"),
        _( "motif_type|m", type="select", options=['all'] + motifs_dict.keys(),
            default="all", help="motif_type: all, one motif or more motifs "
            "(seperated with a ','); available motifs: all, "+
            ', '.join(motifs_dict.keys())+"" ),
        _( "plot", type="select", options=["", "big", "sep"],
            default="", group="plot", help="residue(C3-C'')-plot in one picture "
            "(= big, 6in1), separatly (=sep) or no plot (=default)" ),
        _( "m2f", type="select", options=["True", "False"],
            default="True", help="save the motifs in a txt-file for e.g. superpose"),
        _( "l2f", type="select", options=["True", "False"],
            default="True", help="get the fasta-files from the inputs for multiple alignment"),
        _( "rmsdrange", type="text", default='',
            help="e.g. 1.4,2.3 for upper and under bound")
    ]
    out = [
         _( "jspt_file", file="color_motif.jspt" )
    ]
    RecordsClass = MotifRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "motif.provi"
    def _init( self, *args, **kwargs ):
        self.motif_type = self.motif_type.split(",")
        self.rmsdrange = self.rmsdrange.split(",")
        if self.motif_type[0]=='all':
            self.motif_type=motifs_dict
        self._init_records( utils.path.stem( self.inputs ), **kwargs )
        self._init_parallel( self.inputs, **kwargs )
    def func( self ):
        self.records = find_motifs(self.inputs, self.motif_type, self.m2f, self.l2f, self.rmsdrange ) 
        self.write()
    def _post_exec( self ):
        if self.plot!="":
            plot_found_motifs( self.records, self.plot, self.motif_type, self.output_dir )
        self._make_jspt_color_motif()
        self._make_provi_file(
            pdb_file=self.relpath( self.inputs ),
            jspt_file=self.relpath( self.jspt_file )
        )
        if self.l2f=='True':
            with open(self.output_dir+strftime("%Y-%m-%d_%H-%M-%S", gmtime())+'_motif_fasta.fa', 'w') as fi:
                records_dict = _motif_out_record( self.records )
                for i, motif in enumerate(records_dict):
                    for r in records_dict[motif]:
                        fi.write(r.l2f_str)
        if self.m2f=='True':
            with open(self.output_dir+strftime("%Y-%m-%d_%H-%M-%S", gmtime())+'_all_motifs.txt', 'w') as fi:
                records_dict = _motif_out_record( self.records )
                for i, motif in enumerate(records_dict):
                    for r in records_dict[motif]:
                        fi.write(r.m2f_str)
    def _make_jspt_color_motif( self ):
        records_dict = _motif_out_record( self.records )
        with open( self.jspt_file, "w" ) as fp:
            n = len(records_dict)
            for i, motif in enumerate(records_dict):
                residues = []
                for r in records_dict[motif]:
                    residues.append( "%s:%s" % (r.resno, r.chain) )
                s = "select {%s}; color @{color('roygb', 0, %s, %s)}; wireframe 0.3; \n" % (
                    " or ".join( residues ), n-1, i
                )
                fname2 = self.outpath( "%s_%s.jspt" % (
                    utils.path.stem( self.jspt_file ), motif
                ) )
                with open( fname2 , "w" ) as fp2:
                    fp2.write(
                        "select *; color cpk; wireframe off; " + s +
                        "select none; background black;"
                    )
                fp.write( s )
            fp.write( "select none; background black;" )
       
    


class Superpose( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """A tool to superpose all found motifs. """
    args = [
        _( "motifs_file", type="text" ,
            help="input: txt-file e.g. '~/Downloads/all-motifs.txt'"),
        _( "pdb_repos", type="text",
            help="input: directory, e.g. '~/Downloads'"),
        _( "ava", type="text", default="True",
            help="al vs all (=True), or 1 vs rest "
            "(e.g. ='~/Downloads/3DQB.pdb,A,23'=pdb-file,chain,residue of Ccap)"),
        _( "rmsdrange", type="text", default='0,0',
            help="e.g. 1.4,2.3 for upper and under bound")
    ]
    out = [
        _( "jspt_file", file="color_motif.jspt" )
    ]
    RecordsClass = SuperposeRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "motif.provi"
    def _init( self, *args, **kwargs ):
        self.rmsdrange = map(float, self.rmsdrange.split(","))
        self._init_records( utils.path.stem( self.motifs_file ), **kwargs )
        self._init_parallel( self.motifs_file, **kwargs )
    def func( self ):
        if self.ava=="True":
            self.records = supallall(self.motifs_file, self.pdb_repos, self.rmsdrange, self.output_dir )
        else:
            pass
            #self.records = superposeall()
        #self.write()


    
class StructureGetter( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """A tool to get structure informations. """
    args = [
        _( "inputs", type="text" ,
            help="input: file, files (seperated with a ',') or one directory, "
            "e.g. '~/Downloads' or '~/Downloads/3DQB.pdb,~/Downloads/3SN6'"),
        _( "residue", type="text", help="chain:atom"),
        _( "plot", type="select", options=["", "big", "sep"],
            default="", group="plot", help="residue(C3-C'')-plot in one picture "
            "(= big, 6in1), separatly (=sep) or no plot (=default)" )
    ]
    out = [
        _( "jspt_file", file="color_motif.jspt" )
    ]
    RecordsClass = StrucRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "motif.provi"
    def _init( self, *args, **kwargs ):
        self.residue = self.residue.split(",")
        self._init_records( utils.path.stem( self.inputs ), **kwargs )
        self._init_parallel( self.inputs, **kwargs )
    def func( self ):
        self.records = get_info(self.inputs, self.residue )
        self.write()
    def _post_exec( self ):
        self._make_jspt_color_motif()
        self._make_provi_file(
            pdb_file=self.relpath( self.inputs ),
            jspt_file=self.relpath( self.jspt_file )
        )
        records_list = map( 
            operator.methodcaller( "_asdict" ), 
            self.records 
        )
        print json.dumps( records_list, indent=4 )
        self._plot_all_json()
    def _plot_all_json(self ):
        if self.plot!="":
            ploti(self.output_dir, 'alphaL', self.plot)
    def _make_jspt_color_motif( self ):
        records_dict = _motif_out_struc_record( self.records )
        with open( self.jspt_file, "w" ) as fp:
            n = len(records_dict)
            for i, resnum in enumerate(records_dict):
                residues = []
                for r in records_dict[resnum]:
                    residues.append( "%s:%s" % (r.resno, r.chain) )
                s = "select {%s}; color @{color('roygb', 0, %s, %s)}; wireframe 0.3; \n" % (
                    " or ".join( residues ), n-1, i
                )
                fname2 = self.outpath( "%s_%s.jspt" % (
                    utils.path.stem( self.jspt_file ), resnum
                ) )
                with open( fname2 , "w" ) as fp2:
                    fp2.write(
                        "select *; color cpk; wireframe off; " + s +
                        "select none; background black;"
                    )
                fp.write( s )
            fp.write( "select none; background black;" )