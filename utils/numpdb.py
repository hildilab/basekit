from __future__ import division

import os
import operator
import functools
import numpy as np

from utils import try_int, get_index
from math import vec_dihedral, vec_mag


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
        if self._test_line( line ):
            self._list.append( self._parse_line( line ) )
    def _test_line( self, line ):
        return len( line ) > 0
    def _parse_line( self, line ):
        return line
    def get( self ):
        return np.array( self._list )

class SimpleAtomParser( SimpleParser ):
    def _test_line( self, line ):
        return line.startswith( "ATOM" ) and line[16] in [" ", "A"]

class CoordsParser( SimpleAtomParser ):
    def _parse_line( self, line ):
        return map( float, [ line[30:38], line[38:46], line[46:54] ] )

class ChainParser( SimpleAtomParser ):
    def _parse_line( self, line ):
        return line[21]

class ResnoParser( SimpleAtomParser ):
    def _parse_line( self, line ):
        return int( line[22:26] )

class ResnameParser( SimpleAtomParser ):
    def _parse_line( self, line ):
        return line[17:20]

class AtomnameParser( SimpleAtomParser ):
    def _parse_line( self, line ):
        return line[12:16]

class AltlocParser( SimpleAtomParser ):
    def _parse_line( self, line ):
        return line[16]

class SstrucParser( SimpleParser ):
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
            "coords": CoordsParser(),
            "chain": ChainParser(),
            "resno": ResnoParser(),
            "resname": ResnameParser(),
            "atomname": AtomnameParser(),
            "altloc": AltlocParser(),
            "sstruc": SstrucParser()
        }
    def _parse( self ):
        with open( self.pdb_path, "r" ) as fp:
            for line in fp:
                for p in self._parsers.values():
                    p( line )
        for name in self._parsers.keys():
            p=self._parsers.get( name )
            if p:
                self.__dict__[ "_"+name ] = p.get()
        self.length = len( self._coords )
    def _sele( self, chain=None, resno=None, atomname=None ):
        sele = np.ones( self.length, bool )
        if chain:
            sele &= self._chain==chain
        if resno!=None:
            if isinstance( resno, (list, tuple) ):
                sele &= (self._resno>=resno[0]) & (self._resno<=resno[1])
            else:
                sele &= self._resno==resno
        if atomname:
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
    def _idx( self, **sele ):
        return np.nonzero( self._sele( **sele ) )[0]
    def _idx_first( self, **sele ):
        return self._idx( **sele )[0]
    def iter_chain( self, **sele ):
        for chain in np.unique( self._chain[ self._sele( **sele ) ] ):
            yield chain
    def iter_resno( self, **sele ):
        for chain in self.iter_chain( **sele ):
            sele["chain"] = chain
            sele["resno"] = "CA"
            for resno, chain in np.unique( self._resno[ self._sele( **sele ) ] ):
                print resno
                yield resno, chain
    def __calc_phi_psi( self ):
        self._phi = np.ones( self.length, float )
        self._psi = np.ones( self.length, float )
        for resno, chain in self.iter_resno():
            try:
                idx_ca = self._idx_first( chain=chain, resno=resno, atomname="CA" )
                idx_n = self._idx_first( chain=chain, resno=resno, atomname="N" )
                idx_c = self._idx_first( chain=chain, resno=resno, atomname="C" )
                try:
                    idx_c1 = self._idx_first( chain=chain, resno=resno-1, atomname="C" )
                    phi = vec_dihedral(
                        self._coords[ idx_c1 ], self._coords[ idx_n ],
                        self._coords[ idx_ca ], self._coords[ idx_c ]
                    )
                except:
                    phi = np.nan
                try:
                    idx_n1 = self._idx_first( chain=chain, resno=resno+1, atomname="N" )
                    angle = vec_dihedral(
                        self._coords[ idx_n ], self._coords[ idx_ca ],
                        self._coords[ idx_c ], self._coords[ idx_n1 ]
                    )
                except:
                    psi = np.nan
            except:
                phi = np.nan
                psi = np.nan
            self._phi[ self._sele( chain=chain, resno=resno ) ] = phi
            self._psi[ self._sele( chain=chain, resno=resno ) ] = psi
    def phi( self, **sele ):
        if not hasattr( self, "_phi" ): self.__calc_phi_psi()
        return self._phi[ self._sele( **sele ) ]
    def psi( self, **sele ):
        if not hasattr( self, "_psi" ): self.__calc_phi_psi()
        return self._psi[ self._sele( **sele ) ]
    def dist( self, sele1, sele2 ):
        c1 = self.coords( **sele1 )
        c2 = self.coords( **sele2 )
        v1 = np.sum( c1, axis=0 ) / len(c1)
        v2 = np.sum( c2, axis=0 ) / len(c2)
        return vec_mag( v1 - v2 )



