from __future__ import division

import os
import operator
import functools
import itertools
import collections
from collections import defaultdict
import numpy as np
from cStringIO import StringIO


from utils import try_int, get_index
from math import dihedral, vec_dihedral, mag


# https://pypi.python.org/pypi/Bottleneck
# http://stutzbachenterprises.com/blist/


# http://sourceforge.net/p/pymmlib/code/HEAD/tree/trunk/pymmlib/mmLib/Superposition.py
# https://github.com/bryan-lunt/pdb_tools
# https://github.com/biopython/biopython/tree/master/Bio/PDB


# http://docs.scipy.org/doc/numpy/user/basics.io.genfromtxt.html
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html
# TODO: using genfromtxt would require some preprocessing i.e. all non atom lines
#   need to be removed to parse the atom records; parsing a different
#   record then reqires another call to genfromtxt
# IDEA: option to read only backbone/ mainchain atoms


# PDB FORMAT SPECIFICATIONS
# http://deposit.rcsb.org/adit/docs/pdb_atom_format.html#ATOM
# http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html

HELIX = 1
SHEET = 2


AA1 = {
    'HIS': 'H',
    'ARG': 'R',
    'LYS': 'K',
    'ILE': 'I',
    'PHE': 'F',
    'LEU': 'L',
    'TRP': 'W',
    'ALA': 'A',
    'MET': 'M',
    'PRO': 'P',
    'CYS': 'C',
    'ASN': 'N',
    'VAL': 'V',
    'GLY': 'G',
    'SER': 'S',
    'GLN': 'Q',
    'TYR': 'Y',
    'ASP': 'D',
    'GLU': 'E',
    'THR': 'T'
}
AA3 = dict((v,k) for k, v in AA1.iteritems())



class SimpleParser():
    def __init__( self ):
        self._list = []
    def __call__( self, line ):
        self._list.append( self._parse_line( line ) )
    def _test_line( self, line ):
        return len( line ) > 0
    def _parse_line( self, line ):
        return line
    def get( self ):
        if hasattr( self, "type" ):
            return np.array( self._list, dtype=self.type )
        return np.array( self._list )



pdb_delimiter=(6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6)
pdb_dtype=[
    ('record', '|S6'),          # 0
    ('atomno', np.int),         # 1
    ('empty1', '|S1'),          # 2
    ('atomname', '|S4'),        # 3
    ('altloc', '|S1'),          # 4
    ('resname', '|S3'),         # 5
    ('empty2', '|S1'),          # 6
    ('chain', '|S1'),           # 7
    ('resno', np.int),          # 8
    ('insertion', '|S1'),       # 9
    ('empty3', '|S3'),          # 10
    ('x', np.float),            # 11
    ('y', np.float),            # 12
    ('z', np.float),            # 13
    ('occupancy', np.float),    # 14
    ('bfac', np.float)          # 15
]
# pdb_usecols=(3,4,5,7,8,11,12,13)
pdb_usecols=(3,4,7,8,11,12,13)
pdb_cols=( 
    (0,6), (6,11), (11,12), (12,16), (16,17), (17,20), (20,21), (21,22), 
    (22,26), (26,27), (27,30), (30,38), (38,46), (46,54), (54,60), (60,66)
)



class SstrucParser( SimpleParser ):
    def __init__( self ):
        self._list = []
    def __call__( self, line ):
        self._list.append( self._parse_line( line ) )
    def _test_line( self, line ):
        return line.startswith("HELIX") or line.startswith("SHEET")
    def _parse_line( self, line ):
            if line.startswith("HELIX"):
                return [
                    HELIX,
                    line[19],                   # chain 1
                    try_int( line[21:25] ),     # resno 1
                    line[31],                   # chain 2
                    try_int( line[33:37] ),     # resno 2
                    try_int( line[38:40] ),     # subtype
                    None                        # padding...
                ]
            elif line.startswith("SHEET"):
                return [
                    SHEET,
                    line[21],                   # chain 1
                    try_int( line[22:26] ),     # resno 1
                    line[32],                   # chain 2
                    try_int( line[33:37] ),     # resno 2
                    try_int( line[38:40] ),     # strand sense (subtype)
                    try_int( line[65:69] ),     # resno hbond prev strand
                ]
    def get( self ):
        self._list.sort( key=operator.itemgetter(1,2) )
        return self._list



def sstruc2jmol( sstruc ):
    ret = ""
    for i, ss in enumerate( sstruc ):
        ret += "draw ID 'v%i' vector {%s} {%s};" % (
            i,
            "%0.2f %0.2f %0.2f" % tuple(ss[8]),
            "%0.2f %0.2f %0.2f" % tuple(ss[7])
        )
    return ret

def lsq(y):
    # y = mx + c
    x = np.arange( len(y) )
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq( A, y )[0]
    return [ m*x[0]+c, m*x[-1]+c ]

def axis( coords ):
    return np.array([ lsq(coords[...,i]) for i in range(3) ]).T



class NumPdb:
    atomname_dict = { 
        "CA": " CA ",
        "N": " N  ",
        "C": " C  ",
        "O": " O  ",
        "backbone": ( " CA ", " N  ", " C  ", " O  " ),
        "mainchain": ( " CA ", " N  ", " C  " )
    }
    def __init__( self, pdb_path, features=None ):
        self.pdb_path = pdb_path
        self.features = features
        self._set_parsers()
        self._parse()
    def _set_parsers( self ):
        self._parsers = {
            #"sstruc": SstrucParser()
        }
    def _parse( self ):
        cols = []
        types = []
        for j, c in enumerate(pdb_usecols):
            cols.append( pdb_cols[c] )
            types.append( pdb_dtype[c] )

        types += [ ( 'phi', np.float ), ( 'psi', np.float ) ]

        atoms = []
        header = []

        with open( self.pdb_path, "r" ) as fp:
            for line in fp:
                if line.startswith("ATOM"):# or line.startswith("HETATM"):
                    atoms.append( tuple([ line[ c[0]:c[1] ] for c in cols ] + [ np.nan, np.nan ] ) )
                else:
                    header.append( line )
        
        self._atoms = np.array(atoms, dtype=types)

        # for name, dtype in types:
        #     self.__dict__[ "_"+name ] = atoms[ name ]
        self._coords = np.vstack(( self._atoms['x'], self._atoms['y'], self._atoms['z'] )).T
        # print np.may_share_memory( self._atoms, self._coords )
        # print self._atoms.flags.owndata, self._coords.flags.owndata
        self.length = len( atoms )
        self.__calc_phi_psi()
        self.__tmp_sele = np.ones( self.length, bool )
    def _sele( self, chain=None, resno=None, atomname=None, pre_sele=None, copy=False ):
        if copy:
            sele = np.ones( self.length, bool )
        else:
            self.__tmp_sele.fill( True )
            sele = self.__tmp_sele
        if pre_sele!=None:
            sele &= pre_sele
        if chain!=None:
            if isinstance( chain, (list, tuple) ):
                sele &= reduce( 
                    lambda x, y: x | (self._chain==y), 
                    chain[1:],
                    self._chain==chain[0]
                )
            else:
                sele &= self._chain==chain
        if resno!=None:
            if isinstance( resno, (list, tuple) ):
                sele &= (self._resno>=resno[0]) & (self._resno<=resno[1])
            else:
                sele &= self._resno==resno
        if atomname!=None:
            atomname = self.atomname_dict.get( atomname, atomname )
            if isinstance( atomname, (list, tuple) ):
                sele &= reduce( 
                    lambda x, y: x | (self._atomname==y), 
                    atomname[1:],
                    self._atomname==atomname[0]
                )
            else:
                sele &= self._atomname==atomname
        return sele
    def coords( self, **sele ):
        return self._coords[ self._sele( **sele ) ]
    @property
    def sstruc( self ):
        for i, ss in enumerate( self._sstruc ):
            a = self.axis( sstruc=i )
            ss.extend([ a[1]-a[0], a[0], a[1], i+1 ])
        return self._sstruc
    def axis( self, sstruc=None, **sele ):
        ss = get_index( self._sstruc, sstruc )
        if ss:
            sele["chain"] = ss[1]
            sele["atomname"] = "CA"
            sele["resno"] = [ ss[2], ss[4] ]
        return axis( self.coords( **sele ) )
    def sequence( self, **sele ):
        return "".join([ AA1.get( r, "?" ) for r in self._resname[ self._sele( **sele ) ] ])
    def iter_chain( self, **sele ):
        chain = self._atoms['chain'][0]
        k = 0
        l = 0
        for a in self._atoms:
            if chain!=a['chain']:
                yield self._atoms[k:l], self._coords[k:l]
                chain = a['chain']
                k = l
            l += 1
        yield self._atoms[k:l], self._coords[k:l]
    def iter_resno( self ):
        for atoms, coords in self.iter_chain():
            resno = atoms['resno'][0]
            k = 0
            l = 0
            begin = True
            for a in atoms:
                if resno!=a['resno']:
                    yield atoms[k:l], coords[k:l], begin
                    k = l
                    # detect chain breaks
                    begin = True if (resno!=a['resno']-1) else False
                    resno = a['resno']
                l += 1
            yield atoms[k:l], coords[k:l], begin
    def iter_resno2( self, window ):
        # TODO assumes the first a has a[2]=True and
        #   the first results window has no additional a[2]=True
        it = self.iter_resno()
        for a in it:
            if a[2]:
                result = (a,) + tuple(itertools.islice(it, window-1))
            else:
                result = result[1:] + (a,)
            yield result
    def __calc_phi_psi( self ):

        def get_coords( atoms, coords ):
            return [ coords[ atoms['atomname']==atomname ][0] for atomname in [ ' N  ', ' CA ', ' C  ' ] ]

        for resno3 in self.iter_resno2( 3 ):
            atoms_prev, coords_prev, begin_prev = resno3[0]
            atoms, coords, begin = resno3[1]
            atoms_next, coords_next, begin_next = resno3[2]

            if begin_prev:
                coords_n_prev, coords_ca_prev, coords_c_prev = get_coords( atoms_prev, coords_prev )
                coords_n, coords_ca, coords_c = get_coords( atoms, coords )
                coords_n_next, coords_ca_next, coords_c_next = get_coords( atoms_next, coords_next )
                atoms_prev['psi'] = dihedral( coords_n_prev, coords_ca_prev, coords_c_prev, coords_n )
                atoms['phi'] = dihedral( coords_c_prev, coords_n, coords_ca, coords_c )
            else:
                coords_n, coords_ca, coords_c = coords_n_ca_c
                coords_n_next, coords_ca_next, coords_c_next = get_coords( atoms_next, coords_next )
            
            atoms_next['phi'] = dihedral( coords_c, coords_n_next, coords_ca_next, coords_c_next )
            atoms['psi'] = dihedral( coords_n, coords_ca, coords_c, coords_n_next )
            coords_n_ca_c = ( coords_n_next, coords_ca_next, coords_c_next )
        #for a in self._atoms: print a
    def dist( self, coords1, coords2 ):
        v1 = np.sum( coords1, axis=0 ) / len(coords1)
        v2 = np.sum( coords2, axis=0 ) / len(coords2)
        return mag( v1 - v2 )





