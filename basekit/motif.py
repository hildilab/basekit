from __future__ import with_statement
from __future__ import division

from utils.tool import PyTool


import numpy as np
from numpy import array
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.pyplot import figure
import string, os, sys, argparse
from time import gmtime, strftime
import basekit.utils.numpdb as numpdb
from basekit.utils.timer import Timer
import functools

class AlphaBetaTest():
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
        R Aurora, GD Rose;
            Helix capping (1998), Protein Sci. 7(1):21-38
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
    "alpha_beta": AlphaBetaTest(),
    "alpha_L": AlphaLTest(),
    "alpha_L-more_loose": AlphaLTest_moreloose()
}

def find_all( tool ):
    print tool.output_dir, tool.pdb_file
    print tool.motif_type

def plot_all( inputdir, output_dir, residueplot_in_one, tmp_clearer, motif_type ):
    pdblist=[]
    for dir, subdir, files in os.walk(inputdir):
        for filee in files:
            if filee[-3:]=='pdb':
                pdblist.append(os.path.join(dir,filee))
    plots_all_motifs(
        pdblist, outdir=output_dir, residueplot_in_one=residueplot_in_one,
        tmp_clearer=tmp_clearer, premotif=motif_type
    )

def plots_all_motifs(pdbs, outdir, residueplot_in_one=True, tmp_clearer=True, premotif='all'):
    with Timer("get plots for all motifs"):
        with Timer("file_savings"):
            for fi in pdbs:
                try:
                    print fi
                    npdb=get_npdb(fi)
                except:
                    print 'Error:', fi 
                    pass
                save_found_motifs(fi,npdb,outdir)
        if premotif[0]=='all':
            for motif in motifs_dict:
                with Timer("get plots for motif "+ motif):
                    plot_found_motifs(pdbs, motif, outdir, residueplot_in_one)
        else:
            for motif in premotif:
                with Timer("get plots for motif "+ motif):
                    plot_found_motifs(pdbs, motif, outdir, residueplot_in_one)
    if tmp_clearer: clear(outdir)
    pass

def get_npdb(infile):
    return numpdb.NumPdb( infile)

def save_found_motifs(infile, npdb, outdir):
    name =  os.path.split(infile)[-1]
    tmpdir=outdir+'temp_motifsearch/'
    try:
        os.mkdir(tmpdir)
    except:
        pass
    for motif in motifs_dict:
        try:
            motifs=find_all_motifs(npdb)[motif]
            if motifs!=[]:
                with open(tmpdir+name+motif+'_motif.tmp','w') as fp:
                    for elem in motifs:
                        d = elem[0:2] + tuple(elem[2][0])  + tuple(elem[2][1])
                        fp.write( "%s\n" % ( " ".join( map( str, d ) ) ) )
        except:
            print 'motif_Error', infile
            pass
    pass

def find_all_motifs( npdb ):
    motif_data = {}
    for name, motif_obj in motifs_dict.iteritems():
        motif_data[ name ] = []
        for numatoms in npdb.iter_resno2( motif_obj.length ):
            if motif_obj.test(*numatoms):
                motif_data[ name ].append( motif_obj.make_list(*numatoms) )
    return motif_data

def plot_found_motifs( files, motif, outdir, plot_in_one ):
    n = 0; m = 0; axes=None
    if plot_in_one: fig, axes  = plt.subplots( nrows=2, ncols=3, figsize=(15,15) )
    for index, pos in enumerate(['C3','C2','C1','CCap','C\'','C\'\'']):
        phi=[]; psi=[];
        for p in files:
            try:
                motifs=get_found_motifs(p,motif,outdir)
                for elem in motifs:
                    phi=np.hstack((phi,elem[2][0][index]))
                    psi=np.hstack((psi,elem[2][1][index]))  
            except:
                pass
        if not (phi==[] or psi==[]):
            sele=(np.isnan(phi)==False) & (np.isnan(psi)==False)
            phi = phi[sele]
            psi = psi[sele]
            rama_plot( phi, psi, pos,motif,plot_in_one,n, m,axes, outdir )
            m=m+1
            if m==3:
                m=0
                n=n+1
    if plot_in_one:
        fig.tight_layout()
        name=outdir+motif+'-motif_all'+'_'+strftime("%Y-%m-%d_%H-%M-%S", gmtime())
        plt.suptitle('Ramachandrian plot for '+motif+' motif')
        plt.savefig(name, dpi=300)
        plt.close()
    pass
   
def get_found_motifs(infile,motif,outdir):
    motifs=[]
    tmpdir=os.path.join(outdir,'temp_motifsearch/')
    name =  os.path.split(infile)[-1]
    with open(tmpdir+name+motif+'_motif.tmp','r') as fp:
        for line in fp:
            phi=[];psi=[]
            for pos in range(0,6):
                phi=np.hstack((phi,float(line.split()[pos+2])))
                psi=np.hstack((psi,float(line.split()[pos+8])))
            motifs.append((int(line.split()[0]),line.split()[1],(phi,psi)))
    return motifs

def rama_plot( phi, psi, pos, motif, plot_in_one, n, m, axes, *args ):
    hist2d, phiedges, psiedges = np.histogram2d( phi, psi, bins=180, range=[[-180,180],[-180,180]] )
    extent=[-180,180,-180,180]
    if plot_in_one:
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

def clear(outdir,):
    tmpdir=os.path.join(outdir,'temp_motifsearch/')
    for dir, subdir, files in os.walk(tmpdir):
        for file in files:
            os.remove(dir+file)
    os.rmdir(tmpdir)


class CapsMotifFinder( PyTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "motif_type", "type": "text", "default_value": "alpha_beta" }
    ]
    def _init( self, inputdir, motif_type="alpha_beta", **kwargs ):
        self.inputdir = os.path.abspath( inputdir )
        
        self.motif_type = motif_type.split(",")
        self.func = find_all
        self.output_files = [ "foo.txt" ]


class CapsMotifPlotter( PyTool ):
    args = [
        { "name": "inputdir", "type": "text" },
        { "name": "residueplot_in_one", "type": "checkbox", "default_value": True },
        { "name": "tmp_clearer", "type": "checkbox", "default_value": True },
        { "name": "motif_type", "type": "select", "options": ['all'] + motifs_dict.keys(), "default_value": "all" }
    ]
    def _init( self, inputdir, residueplot_in_one=True, tmp_clearer=True, motif_type="all", **kwargs ):
        self.inputdir = os.path.abspath( inputdir )
        self.residueplot_in_one=residueplot_in_one
        self.tmp_clearer=tmp_clearer
        self.motif_type = motif_type.split(",")
        print self.motif_type
        self.output_files = [ "foo.txt" ]
    def func( self ):
        plot_all( self.inputdir, self.output_dir, self.residueplot_in_one, self.tmp_clearer, self.motif_type )
    

