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






DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "motif" )



MotifRecord = collections.namedtuple( 'MotifRecord', [
    'motif', 'pdb_id', 'chain', 'resno', 'rmsd',
    'phi_angles', 'psi_angles', 'no'
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








motifs_dict = {
    "alpha_beta_motif": AlphaBetaMotifTest(),
    "alpha_L_motif": AlphaLTest(),
    "alpha_L_motif-fitted": AlphaLTest_fitted(),
    "st_motif": STTest()
}


def _superpose_test( numa, sele_motif ):

    numa_ref=numpdb.NumPdb("/home/student/Johanna/Projekte/caps/cap-referenceset/2REB.pdb")
    chain=sele_motif['chain']
    cd=sele_motif['cd']
    sele_ref = {
            "chain": "A", "resno": [ 86-4, 86 ]
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



def supallall( motifs_file, pdb_repos, *args ):
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
    hist2d, phiedges, psiedges = np.histogram2d( arr, arr, bins=25,
                                                range=[[0,0.25],[0,0.25]]
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



def _make_record( motif, name, motif_line, i, rmsd ):
    return MotifRecord(
        motif, name, motif_line[1], motif_line[0], rmsd,
        list(motif_line[2][0]), list(motif_line[2][1]), i
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

def find_motifs( infile, motif_type, m2f, l2f ):
    # generate a numpy-pdb
    # get the motif-dicc and parses the motifs to Motif-Record
    records = []
    if l2f: fi=open('/home/student/Johanna/Projekte/caps/fasta.fa', 'a')
    if m2f: fo=open('/home/student/Johanna/Projekte/caps/motifs.txt', 'a')
    npdb=numpdb.NumPdb(infile)
    pdb_id = utils.path.stem( infile )
    for motif in motif_type:
        motif_found=_find_all_motifs(npdb)[motif]
        if motif_found!=[]:
            for i, elem in enumerate(motif_found):
                rmsd=_superpose_test( npdb, {'chain':elem[1], 'cd':elem[0]} )
                resstr=""
                if l2f:
                    for i in range(elem[0]-4, elem[0]):
                        resstr=resstr+bio.AA1.get(npdb.iter_resno2(2)[i-4][0][0]["resname"], "X")
                    fi.write(">"+pdb_id+"|"+elem[1]+str(elem[0])+"\n")
                    fi.write(resstr)
                    fi.write('\n')
                if m2f:
                    fo.write(pdb_id+" "+elem[1]+" "+str(elem[0])+"\n")
                #if rmsd > 0.3 and rmsd < 0.4:
                #    print pdb_id, rmsd, elem[1], elem[0]
                #    records.append(_make_record(motif, pdb_id, elem, i, rmsd))
                #if rmsd < 0.4:
                #    records.append(_make_record(motif, pdb_id, elem, i, rmsd))
                records.append(_make_record(motif, pdb_id, elem, i, rmsd))

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
            default="True", help="save the motifs in a txt-file"),
        _( "l2f", type="select", options=["True", "False"],
            default="True", help="get the fasta-files from the inputs")
    ]
    out = [
         _( "jspt_file", file="color_motif.jspt" )
    ]
    RecordsClass = MotifRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "motif.provi"
    def _init( self, *args, **kwargs ):
        self.motif_type = self.motif_type.split(",")
        if self.motif_type[0]=='all':
            self.motif_type=motifs_dict
        self._init_records( utils.path.stem( self.inputs ), **kwargs )
        self._init_parallel( self.inputs, **kwargs )
    def func( self ):
        self.records = find_motifs(self.inputs, self.motif_type, self.m2f, self.l2f ) 
        self.write()
    def _post_exec( self ):
        if self.plot!="":
            plot_found_motifs( self.records, self.plot, self.motif_type, self.output_dir )
        self._make_jspt_color_motif()
        self._make_provi_file(
            pdb_file=self.relpath( self.inputs ),
            jspt_file=self.relpath( self.jspt_file )
        )
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
            "(e.g. ='~/Downloads/3DQB.pdb,A,23'=pdb-file,chain,residue of Ccap)")
    ]
    out = [
        _( "jspt_file", file="color_motif.jspt" )
    ]
    RecordsClass = SuperposeRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "motif.provi"
    def _init( self, *args, **kwargs ):
        self._init_records( utils.path.stem( self.motifs_file ), **kwargs )
        self._init_parallel( self.motifs_file, **kwargs )
    def func( self ):
        if self.ava=="True":
            self.records = supallall(self.motifs_file, self.pdb_repos, self.output_dir )
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