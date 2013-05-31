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
from utils import dir_walker, boolean
import operator


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



def find_all( inputdir, output_dir, tmp_clearer, motif_type, printer ):
    if os.path.isdir(inputdir[0]):
        pdblist = map( operator.itemgetter(1), dir_walker(inputdir[0],'.+\.pdb') )
    else:
        pdblist = inputdir
    find_and_save_motifs(pdblist, motif_type, output_dir, printer )
    if tmp_clearer: clear(output_dir)
    pass

def plot_all( inputdir, output_dir, residueplot_in_one, tmp_clearer, motif_type, printer ):
    motif_ft={}
    if os.path.isdir(inputdir[0]):
        pdblist = map( operator.itemgetter(1), dir_walker(inputdir[0],'.+\.pdb') )
    else:
        pdblist = inputdir
    with Timer("get plots for all motifs"):
        motif_ft=find_and_save_motifs(pdblist, motif_type, output_dir, printer )
        if motif_type[0]=='all':
            for motif in motifs_dict:
                if motif_ft[motif]:
                    with Timer("get plots for motif "+ motif):
                        plot_found_motifs(pdblist, motif, output_dir, residueplot_in_one)
                else: print '### No '+motif+'-motif found!'
        else:
            for motif in motif_type:
                if motif_ft[motif]:
                    with Timer("get plots for motif "+ motif):
                        plot_found_motifs(pdblist, motif, output_dir, residueplot_in_one)
                else: print '### No '+motif+'-motif found!'
    if tmp_clearer: clear(output_dir)
    pass

def find_and_save_motifs( pdblist, motif_type, output_dir, printer ):
    motif_ft={}
    with Timer("find motifs"):
        for fi in pdblist:
            try:
                npdb=numpdb.NumPdb(fi)
                if printer: print "motifs in file "+ fi
                if motif_type[0]=='all': motif_ft = motif_sap(npdb, motifs_dict, printer, output_dir, fi)
                else: motif_ft = motif_sap(npdb, motif_type, printer, output_dir, fi)
            except:
                print 'Error:', fi 
                pass
    return motif_ft

def motif_sap(npdb, motifs, printer, output_dir, fi):
    motif_ft={}
    for motif in motifs:
        motif_list=find_all_motifs(npdb)[motif]
        if motif_list!=[]:
            if printer: print '### '+motif+'-motif: Ccap at residue: ',map(operator.itemgetter(0), motif_list)
            motif_saver(fi, output_dir, motif, motif_list)
            motif_ft[motif]=True
        else:
            if printer: print '###','No '+motif+'-motif found!'
            motif_ft[motif]=False
    return motif_ft


def find_all_motifs( npdb ):
    motif_data = {}
    for name, motif_obj in motifs_dict.iteritems():
        motif_data[ name ] = []
        for numatoms in npdb.iter_resno2( motif_obj.length ):
            if motif_obj.test(*numatoms):
                motif_data[ name ].append( motif_obj.make_list(*numatoms) )
    return motif_data

def motif_saver(infile, outdir, motif, motif_list):
    name =  os.path.split(infile)[-1]
    tmpdir=outdir+'temp_motifsearch/'
    try:
        os.mkdir(tmpdir)
    except:
        pass
    with open(tmpdir+name+motif+'_motif.tmp','w') as fp:
        for elem in motif_list:
            d = elem[0:2] + tuple(elem[2][0])  + tuple(elem[2][1])
            fp.write( "%s\n" % ( " ".join( map( str, d ) ) ) )
    pass

def plot_found_motifs( files, motif, outdir, plot_in_one ):
    n = 0; m = 0; axes=None
    if plot_in_one: fig, axes  = plt.subplots( nrows=2, ncols=3, figsize=(15,15) )
    for index, pos in enumerate(['C3','C2','C1','CCap','C\'','C\'\'']):
        phi=[]; psi=[];
        for p in files:
            try:
                motifs=get_found_motifs(p,motif,outdir)
                if motifs==[]:
                    pass
                else:
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
    """A tool that finds cap motifs. """
    args = [
        { "name": "inputdir", "type": "text" ,
          "help": "input: file, files (seperated with a ',') or one directory, "
                  "e.g. '~/Downloads' or '~/Downloads/3DQB.pdb,~/Downloads/3SN6'"},
        { "name": "tmp_clearer", "type": "checkbox", "default_value": False,
          "help": "deletes the temp-files (=False, default)" },
        { "name": "motif_type", "type": "select", "options": ['all'] + motifs_dict.keys(),
          "default_value": "all", "help": "motif_type: all, one motif or more motifs "
                           "(seperated with a ','); available motifs: all, "+
                           ', '.join(motifs_dict.keys())+"" },
        { "name": "printer", "type": "checkbox", "default_value": True,
          "help": "prints motifs for each file (=True, default)" }
    ]
    def _init( self, inputs, tmp_clearer=False, motif_type="all", printer=True, **kwargs ):
        self.inputs = inputs.split(",")
        self.tmp_clearer=tmp_clearer
        self.motif_type = motif_type.split(",")
        self.printer = printer
        self.output_files = [ "foo.txt" ]
    def func( self ):
        find_all( self.inputs, self.output_dir, self.tmp_clearer, self.motif_type, self.printer )


class CapsMotifPlotter( PyTool ):
    """A tool that finds cap motifs and pots their dihedral angles. """
    args = [
        { "name": "inputs", "type": "text" ,
          "help": "input: file, files (seperated with a ',') or one directory, "
                  "e.g. '~/Downloads' or '~/Downloads/3DQB.pdb,~/Downloads/3SN6'"},
        { "name": "bigplot", "type": "checkbox", "default_value": True,
          "help": "residue(C3-C'')-plot in one picture (=True, default, 6in1) "
                  "or separatly (=False) " },
        { "name": "tmp_clearer", "type": "checkbox", "default_value": True,
          "help": "deletes the temp-files (=True, default)" },
        { "name": "motif_type", "type": "select", "options": ['all'] + motifs_dict.keys(),
          "default_value": "all", "help": "motif_type: all, one motif or more motifs "
                           "(seperated with a ','); available motifs: all, "+
                           ', '.join(motifs_dict.keys())+"" },
        { "name": "printer", "type": "checkbox", "default_value": False,
          "help": "prints motifs for each file (=True, default=False)" }
    ]
    def _init( self, inputs, bigplot=True, tmp_clearer=True, motif_type="all", printer=False, **kwargs ):
        self.inputs = inputs.split(",")
        self.bigplot = bigplot
        self.tmp_clearer=tmp_clearer
        self.motif_type = motif_type.split(",")
        self.printer = printer
        self.output_files = [ "foo.txt" ]
    def func( self ):
        plot_all( self.inputs, self.output_dir, self.bigplot, self.tmp_clearer, self.motif_type, self.printer )
    

