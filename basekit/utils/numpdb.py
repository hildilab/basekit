from __future__ import division

import os
import operator
import functools
import itertools
import collections
import logging

import numpy as np
np.seterr( all="raise" )

import basekit.utils.path
from basekit.utils import try_int, get_index
from math import dihedral, vec_dihedral, mag, axis


logging.basicConfig()
LOG = logging.getLogger('numpdb')
# LOG.setLevel( logging.ERROR )
LOG.setLevel( logging.WARNING )
# LOG.setLevel( logging.DEBUG )



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


RESIDUES = {
    'nucleotides': frozenset([ 
        "  C", "  U", "  G", "  A",
        " DC", " DT", " DG", " DA", " DI", " DU"
    ]),
    'aminoacids': frozenset([
        "TPO", "MSE", "HYP", "SMF", "ALC", "KCX"
    ])
}

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

ATOMS = { 
    "CA": " CA ",
    "CD1": " CD1",
    "CD2": " CD2",
    "N": " N  ",
    "C": " C  ",
    "O": " O  ",
    "backbone": frozenset([ " N  ", " CA ", " C  ", " O  " ]),
    "mainchain": frozenset([ " N  ", " CA ", " C  " ])
}


def numsele( string ):
    sele = { "chain": None, "resno": None, "atomname": None }
    atomname = string.split(".")
    if len(atomname)>1 and atomname[1]:
        sele["atomname"] = atomname[1][0:4]
    chain = atomname[0].split(":")
    if len(chain)>1 and chain[1]:
        sele["chain"] = chain[1][0]
    if chain[0]:
        sele["resno"] = int(chain[0])
    return sele


def numdist( numa1, numa2 ):
    return mag( numa1.center() - numa2.center() )



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
        if hasattr( self, "type" ):
            return np.array( self._list, dtype=self.type )
        return np.array( self._list )



PDB_DELIMITER=(6,5,1,4,1,3,1,1,4,1,3,8,8,8,6,6,6,4,2,2)
PDB_DTYPE=[
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
    ('bfac', np.float),         # 15
    ('empty4', '|S6'),          # 16
    ('segment', '|S4'),         # 17
    ('element', '|S2'),         # 18
    ('charge', '|S2'),          # 19
]
PDB_USECOLS=(1,3,4,5,7,8,9,11,12,13,14,15)
PDB_COLS=( 
#   0        1        2        3        4        5        6        7
    ( 0, 6), ( 6,11), (11,12), (12,16), (16,17), (17,20), (20,21), (21,22),
#   8        9        10       11       12       13       14       15
    (22,26), (26,27), (27,30), (30,38), (38,46), (46,54), (54,60), (60,66),
#   16       17       18       19
    (66,72), (72,76), (76,78), (78,80)
)
PDB_DEFAULTS = {
    "record": "ATOM",
    "atomno": 1,
    "atomname": " C  ",
    "altloc": "",
    "resname": "GLY",
    "chain": "",
    "resno": 1,
    "insertion": "",
    "x": 0.0, "y": 0.0, "z": 0.0,
    "occupancy": 1.0,
    "bfac": 0.0,
    "segment": "",
    "element": "",
    "charge": ""
}
PDB_ATOM_TMPL = "{record:<6}{atomno:>5} {atomname:>4}{altloc:>1}{resname:>3} " \
                "{chain:>1}{resno:>4}{insertion:>1}   {x:8.3f}{y:8.3f}{z:8.3f}" \
                "{occupancy:6.2f}{bfac:6.2f}      {segment:<4}{element:>2}{charge:>2}\n"

def numdefaults( natom, defaults ):
    d = {}
    for k, v in defaults.iteritems():
        try:
            d[k] = natom[k]
        except:
            d[k] = v
    return d

def pdb_line( natom ):
    return PDB_ATOM_TMPL.format( 
        **numdefaults( natom, PDB_DEFAULTS )
    )

def xyzr_line( natom ):
    pass



SstrucRecord = collections.namedtuple( 'SstrucRecord', [
    'type', 'subtype',
    'chain1', 'resno1',
    'chain2', 'resno2',
    'hbond'
])

class SstrucParser( SimpleParser ):
    def __init__( self ):
        self._list = []
    def _test_line( self, line ):
        return line.startswith("HELIX") or line.startswith("SHEET")
    def _parse_line( self, line ):
            if line.startswith("HELIX"):
                return SstrucRecord(
                    HELIX,
                    try_int( line[38:40] ),     # subtype
                    line[19],                   # chain 1
                    try_int( line[21:25] ),     # resno 1
                    line[31],                   # chain 2
                    try_int( line[33:37] ),     # resno 2
                    None                        # padding...
                )
            elif line.startswith("SHEET"):
                return SstrucRecord(
                    SHEET,
                    try_int( line[38:40] ),     # strand sense (subtype)
                    line[21],                   # chain 1
                    try_int( line[22:26] ),     # resno 1
                    line[32],                   # chain 2
                    try_int( line[33:37] ),     # resno 2
                    try_int( line[65:69], False ),     # resno hbond prev strand
                )
    def get( self ):
        self._list.sort( key=operator.attrgetter("chain1", "resno1") )
        return self._list






class NumAtoms:
    def __init__( self, atoms, coords, flag=None ):
        self._atoms = atoms
        self._coords = coords
        self.flag = flag
        self.length = len( atoms )
    def __getitem__( self, key ):
        if key=='xyz':
            return self._coords
        else:
            return self._atoms[ key ]
    def __setitem__( self, key, value ):
        if key=='xyz':
            self._coords = value
        else:
            self._atoms[ key ] = value
    def sele( self, chain=None, resno=None, atomname=None, sele=None, invert=None ):
        atoms = self._atoms
        if sele==None:
            sele = np.ones( self.length, bool )
        if chain!=None:
            if isinstance( chain, collections.Iterable ):
                tmps = atoms['chain']==chain[0]
                for an in chain[1:]: tmps |= atoms['chain']==an
                sele &= tmps
            else:
                sele &= atoms['chain']==chain
        if resno!=None:
            if isinstance( resno, collections.Iterable ):
                sele &= (atoms['resno']>=resno[0]) & (atoms['resno']<=resno[1])
            else:
                sele &= atoms['resno']==resno
        if atomname!=None:
            if isinstance( atomname, collections.Iterable ) and not isinstance( atomname, basestring ):
                atomname = [ ATOMS.get( a, a ) for a in atomname ]
                tmps = atoms['atomname']==atomname[0]
                for an in atomname[1:]: tmps |= atoms['atomname']==an
                sele &= tmps
            else:
                atomname = ATOMS.get( atomname, atomname )
                sele &= atoms['atomname']==atomname
        if invert:
            np.logical_not( sele, sele )
        return sele
    def slice( self, begin, end, flag=None ):
        return NumAtoms( self._atoms[begin:end], self._coords[begin:end], flag=flag )
    def copy( self, **sele ):
        _sele = self.sele( **sele )
        return NumAtoms( self._atoms[ _sele ], self._coords[ _sele ] )
    def _select( self, **sele ):
        coords = self._coords
        atoms = self._atoms
        if len(sele):
            _sele = self.sele( **sele )
            coords = self._coords[ _sele ]
            atoms = self._atoms[ _sele ]
        return coords, atoms
    def get( self, key, **sele ):
        coords, atoms = self._select( **sele )
        if key=='xyz':
            return coords
        else:
            return atoms[ key ]
    def index( self, first=False, last=False, **sele ):
        indices = np.nonzero( self.sele( **sele ) )[0]
        if first:
            return indices[0]
        elif last:
            return indices[-1]
        else:
            return indices
    def iter_chain( self, **sele ):
        coords, atoms = self._select( **sele )
        if len(atoms)==0: return
        chain = atoms['chain'][0]
        k = 0
        l = 0
        for a in atoms:
            if chain!=a['chain']:
                yield self.slice( k, l )
                chain = a['chain']
                k = l
            l += 1
        yield self.slice( k, l )
    def iter_sstruc( self, **sele ):
        for numatoms in self.iter_chain( **sele ):
            sstruc = numatoms['sstruc'][0]
            atoms = numatoms._atoms
            k = 0
            l = 0
            for a in atoms:
                if sstruc!=a['sstruc']:
                    yield numatoms.slice( k, l )
                    k = l
                    sstruc = a['sstruc']
                l += 1
            yield numatoms.slice( k, l )
    def _iter_resno( self, **sele ):
        for numatoms in self.iter_chain( **sele ):
            res = ( numatoms['resno'][0], numatoms['insertion'][0] )
            k = 0
            l = 0
            for a in numatoms._atoms:
                if res!=( a['resno'], a['insertion'] ):
                    numa = numatoms.slice( k, l )
                    yield numa
                    k = l
                    res = ( a['resno'], a['insertion'] )
                l += 1
            yield numatoms.slice( k, l )
    def iter_resno( self, **sele ):
        for numatoms in self.iter_chain( **sele ):
            
            def get_slice( k, l, numa_prev ):
                numa = numatoms.slice( k, l, flag=True )
                if numa_prev:
                    if numa_prev['resno'][0]==numa['resno'][0]-1:
                        flag = False
                    elif numdist( numa_prev.copy( atomname="C" ), numa.copy( atomname="N" ) ) > 1.4:
                        flag = True
                    else:
                        flag = False
                    numa.flag = flag
                return numa

            res = ( numatoms['resno'][0], numatoms['insertion'][0] )
            numa_prev = None
            # atoms = numatoms._atoms
            k = 0
            l = 0

            for a in numatoms._atoms:
                if res!=( a['resno'], a['insertion'] ):
                    numa = get_slice( k, l, numa_prev )
                    yield numa
                    k = l
                    res = ( a['resno'], a['insertion'] )
                    numa_prev = numa
                l += 1
            yield get_slice( k, l, numa_prev )
    def iter_resno2( self, window, **sele ):
        # (TODO) assumes the first a has a.flag==True
        it = self.iter_resno( **sele )
        for a in it:
            if a.flag:
                result = [ a ]
                for b in it:
                    if b.flag: 
                        result = [ b ]
                    else: 
                        result.append( b )
                    if len( result ) == window:
                        break
                if len( result ) == window:
                    result = tuple(result)
                else:
                    return
            else:
                result = result[1:] + (a,)
            yield result
    def axis( self, **sele ):
        return axis( self.get( 'xyz', **sele ) )
    def sequence( self, **sele ):
        return "".join([ AA1.get( a['resname'][0], "?" ) for a in self.iter_resno( **sele ) ])
    def center( self, **sele ):
        coords = self.get( 'xyz', **sele )
        return np.sum( coords, axis=0 ) / len(coords)
    def dist( self, sele1, sele2 ):
        return mag( self.center( **sele1 ) - self.center( **sele2 ) )
    def write( self, file_name, **sele ):
        coords, atoms = self._select( **sele )
        with open( file_name, "w" ) as fp:
            for natom in atoms:
                fp.write( pdb_line( natom ) )





class NumPdb:
    def __init__( self, pdb_path, features=None ):
        self.pdb_path = pdb_path
        self.pdb_id = basekit.utils.path.stem( pdb_path )
        self.features = {
            "phi_psi": True,
            "sstruc": True,
            "backbone_only": False,
            "protein_only": True,
            "detect_incomplete": True
        }
        if features: self.features.update( features )
        self._parse()
    def __getattr__(self, attr):
        # delegate access off nonexistant attributes to the NumAtoms instance
        return getattr(self.numatoms, attr)
    def _parse( self ):
        cols = []
        types = []
        for c in PDB_USECOLS:
            cols.append( PDB_COLS[c] )
            types.append( PDB_DTYPE[c] )

        extra = []
        parsers = {}
        if self.features["detect_incomplete"]:
            types += [ ( 'incomplete', np.bool ) ]
            extra += [ False ]
        if self.features["phi_psi"]:
            types += [ ( 'phi', np.float ), ( 'psi', np.float ) ]
            extra += [ np.nan, np.nan ]
        if self.features["sstruc"]:
            types += [ ( 'sstruc', '|S1' ) ]
            extra += [ " " ]
            parsers[ "sstruc" ] = SstrucParser()

        atoms = []
        atoms_append = atoms.append
        header = []
        header_append = header.append

        with open( self.pdb_path, "r" ) as fp:
            backbone = ATOMS['backbone']
            nucleotides = RESIDUES['nucleotides']
            aminoacids = RESIDUES['aminoacids']
            nbo = not self.features["backbone_only"]
            po = self.features["protein_only"]
            altloc = (' ', 'A', '1')
            keys = ('ATOM  ', 'HETATM', 'MODEL ')
            parsrs = parsers.values()
            tupl = tuple

            for line in fp:
                if line[0:6] in keys:
                    break
                header_append( line )
                for p in parsrs:
                    p( line )
            
            for line in itertools.chain( [line], fp ):
                if line[0:4]=="ATOM" or ( line[0:6]=="HETATM" and line[17:20] in aminoacids ):
                    if po and line[17:20] in nucleotides:
                        continue
                    if ( nbo or line[12:16] in backbone ) and line[16] in altloc:
                        atoms_append( tupl( [ line[ c[0]:c[1] ] for c in cols ] + extra ) )
                elif line[0:5]=="MODEL":
                    # TODO add model field
                    pass
                elif line[0:6]=="CONECT":
                    # stop condition
                    break
                    
        for name in parsers.keys():
            p=parsers.get( name )
            if p:
                self.__dict__[ "_%s" % name ] = p.get()

        atoms = np.array(atoms, dtype=types)
        coords = np.vstack(( atoms['x'], atoms['y'], atoms['z'] )).T
        # print np.may_share_memory( self.atoms, self.coords )
        # print self.atoms.flags.owndata, self.coords.flags.owndata

        self._atoms = atoms
        self._coords = coords
        self._header = header
        self.numatoms = NumAtoms( self._atoms, self._coords )

        self.length = len( atoms )
        if self.features["detect_incomplete"]:
            self.__calc_incomplete()
        if self.features["phi_psi"]:
            self.__calc_phi_psi()
        if self.features["sstruc"]:
            self.__calc_sstruc()
    def __calc_incomplete( self ):
        bb_subset = ATOMS['backbone'].issubset
        try:
            for numa in self._iter_resno():
                if not bb_subset( numa["atomname"] ):
                    numa['incomplete'] = True
                    # print numa._atoms
        except Exception as e:
            LOG.error( "[%s] calc incomplete (%s) => %s" % ( self.pdb_id, e, numa['atomno'][0] ) )
        sele = self._atoms["incomplete"]==False
        if sele!=[]:
            self._atoms = self._atoms[ sele ]
            self._coords = self._coords[ sele ]
            self.numatoms = NumAtoms( self._atoms, self._coords )
    def __calc_sstruc( self ):
        for ss in self._sstruc:
            try:
                idx_beg = self.index( chain=ss.chain1, resno=ss.resno1, first=True )
                idx_end = self.index( chain=ss.chain2, resno=ss.resno2, last=True )
                if ss.type==HELIX:
                    self.slice( idx_beg, idx_end )['sstruc'] = "H"
                elif ss.type==SHEET:
                    self.slice( idx_beg, idx_end )['sstruc'] = "E"
            except Exception as e:
                LOG.error( "[%s] calc sstruc (%s) => %s" % ( self.pdb_id, e, ss ) )
    def __calc_phi_psi( self, dihedral=dihedral ):
        mainchain = ATOMS['mainchain']
        for na_curr, na_next in self.iter_resno2( 2 ):
            try:
                xyz_n, xyz_ca, xyz_c = na_curr.get( 'xyz', atomname=mainchain )
                xyz_n_next, xyz_ca_next, xyz_c_next = na_next.get( 'xyz', atomname=mainchain )
                na_curr['psi'] = dihedral( xyz_n, xyz_ca, xyz_c, xyz_n_next )
                na_next['phi'] = dihedral( xyz_c, xyz_n_next, xyz_ca_next, xyz_c_next )
            except Exception as e:
                LOG.error( "[%s] calc phi/psi (%s) => %s, %s" % (
                    self.pdb_id, e, na_curr['atomname'], na_next['atomname']
                ))
        # for a in self.iter_resno(): print a._atoms[0], a.flag





