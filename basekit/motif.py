from __future__ import with_statement
from __future__ import division




import numpy as np
from numpy import array

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.pyplot import figure

import string, os, sys, argparse, itertools

from time import gmtime, strftime
import functools, operator, collections
import basekit.utils.numpdb as numpdb

from collections import defaultdict

import utils.math, utils.path
from utils import (
    try_int, get_index, boolean, iter_window, memoize, 
    copy_dict, dir_walker
)
from utils.timer import Timer
from utils.db import get_pdb_files, create_table
from utils.tool import (
    PyTool, DbTool, RecordsMixin, ParallelMixin, ProviMixin
)
import utils.path

import utils.numpdb as numpdb








DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "motif" )



MotifRecord = collections.namedtuple( 'MotifRecord', [
    'motif', 'pdb_id', 'chain', 'resno',
    'phi_angles', 'psi_angles', 'no'
])




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

class AlphaLTest_moreloose():
    """
    * H-bond between CO of residue C3 and NH of residue Cd
    * phi angles of residues Ccap=[-99,-81], Cd=[79,91], Cdd=[-113,-87], Cddd=[-112,-88]
    * psi angles of residues Ccap=[-16,-2], Cd=[3,21], Cdd=[110,144]
    Reference:
        self
    """
    length = 6
    def test(self, C3, C2, C1, Ccap, Cd, Cdd):#, Cddd):
        C3_dist = C3.copy(atomname="O")
        Cd_dist = Cd.copy(atomname="N")
        distance1 = numpdb.numdist( C3_dist, Cd_dist )
        return distance1<3.9 and 79<Cd['phi'][0]<91 and 3<Cd['psi'][0]<21 
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










motifs_dict = {
    "alpha_beta_motif": AlphaBetaMotifTest(),
    "alpha_L_motif": AlphaLTest(),
    "alpha_L-more_loose_motif": AlphaLTest_moreloose(),
    "st_motif": STTest()
}





def _make_record( motif, name, motif_line, i ):
    return MotifRecord(
        motif, name, motif_line[1], motif_line[0], 
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


def find_motifs( infile, motif_type ):
    # generate a numpy-pdb
    # get the motif-dicc and parses the motifs to Motif-Record
    records = []
    npdb=numpdb.NumPdb(infile)
    pdb_id = utils.path.stem( infile )
    for motif in motif_type:
        motif_found=_find_all_motifs(npdb)[motif]
        if motif_found!=[]:
            for i, elem in enumerate(motif_found):
                records.append(_make_record(motif, pdb_id, elem, i)) 
    return records


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
            phi=np.hstack((phi,motif_pos.phi_angles))
            psi=np.hstack((psi,motif_pos.psi_angles))
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


class CapsMotifFinder( PyTool, RecordsMixin, ParallelMixin, ProviMixin ):
    """A tool that finds cap motifs. """
    args = [
        { "name": "inputs", "type": "text" ,
          "help": "input: file, files (seperated with a ',') or one directory, "
                  "e.g. '~/Downloads' or '~/Downloads/3DQB.pdb,~/Downloads/3SN6'"},
        { "name": "motif_type", "type": "select", "options": ['all'] + motifs_dict.keys(),
          "default_value": "all", "help": "motif_type: all, one motif or more motifs "
                           "(seperated with a ','); available motifs: all, "+
                           ', '.join(motifs_dict.keys())+"" },
        { "name": "plot", "type": "select", "options": ["", "big", "sep"],
          "default_value": "", "group": "plot",
          "help": "residue(C3-C'')-plot in one picture (= big, 6in1), "
                  "separatly (=sep) or no plot (=default)" }
    ]
    RecordsClass = MotifRecord
    tmpl_dir = TMPL_DIR
    provi_tmpl = "motif.provi"
    def _init( self, inputs, motif_type="all", plot="", **kwargs ):
        self.inputs = self.abspath( inputs )
        self.plot=plot
        self.motif_type = motif_type.split(",")
        if self.motif_type[0]=='all':
            self.motif_type=motifs_dict
        self.jspt_file = self.outpath( "color_motif.jspt" )
        self.output_files = []
        self._init_records( utils.path.stem( inputs ), **kwargs )
        self._init_parallel( self.inputs, **kwargs )
    def func( self ):
        if self.parallel:
            self._make_tool_list()
            tool_results = self._func_parallel()
            self.records = list(itertools.chain.from_iterable(
                map( operator.attrgetter( "records" ), tool_results )
            ))
        else:
            self.records = find_motifs(self.inputs, self.motif_type )
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
    

