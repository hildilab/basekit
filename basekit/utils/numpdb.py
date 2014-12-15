from __future__ import division

import os
import operator
import itertools
import collections
import logging
import json

import numpy as np
np.seterr( all="raise" )

import basekit.utils.path
from basekit.utils.align import aligner
from basekit.utils import (
    try_int, try_float, get_index, iter_window,
    memoize_m, listify
)
from math import mag, axis, Superposition
from bio import AA1

try:
    from cgeom import dihedral
except ImportError as e:
    print "cgeom import error", e
    from math import dihedral


logging.basicConfig()
LOG = logging.getLogger('numpdb')
# LOG.setLevel( logging.ERROR )
LOG.setLevel( logging.WARNING )
# LOG.setLevel( logging.DEBUG )


DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
BASEKIT_DIR = os.path.split( PARENT_DIR )[0]
DATA_DIR = os.path.join( BASEKIT_DIR, "data", "pdb" )


# https://pypi.python.org/pypi/Bottleneck
# http://stutzbachenterprises.com/blist/
# http://wiki.python.org/moin/PythonSpeed/PerformanceTips

# http://sourceforge.net/p/pymmlib/code/HEAD/tree/trunk/pymmlib/mmLib/Superposition.py
# http://sourceforge.net/p/pymmlib/code/1197/tree/
# https://github.com/bryan-lunt/pdb_tools
# https://github.com/biopython/biopython/tree/master/Bio/PDB
# http://cci.lbl.gov/~rwgk/tmp/pdb-v3.2/2009_04_14/news.html
# http://sourceforge.net/p/cctbx/code/17756/tree/trunk/iotbx/pdb/hierarchy.py
# http://www.csb.pitt.edu/ProDy/_modules/prody/proteins/pdbfile.html

# http://docs.scipy.org/doc/numpy/user/basics.io.genfromtxt.html
# http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html
# TODO: using genfromtxt would require some preprocessing
#   i.e. all non atom lines
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
    'aminoacids': frozenset([  ])  # read from file, see below
}
try:
    additional_aminoacids = [ 'MPR', 'PSU', 'PYR', '4AF' ]
    with open( os.path.join( DATA_DIR, "aa.json" ), "r" ) as fp:
        RESIDUES['aminoacids'] = frozenset(
            json.load( fp ) + additional_aminoacids
        )
except:
    LOG.warning( 'Aminoacid list could not be loaded.' )


ATOMS = {
    "CA": " CA ",
    "NZ": " NZ ",
    "NE": " NE ",
    "NE1": " NE1",
    "NE2": " NE2",
    "NH1": " NH1",
    "NH2": " NH2",
    "ND2": " ND2",
    "ND1": " ND1",
    "CZ": " CZ ",
    "CZ2": " CZ2",
    "CZ3": " CZ3",
    "CE1": " CE1",
    "CE2": " CE2",
    "CE3": " CE3",
    "CE": " CE ",
    "CH2": " CH2",
    "OH": " OH ",
    "OE1": " OE1",
    "OE2": " OE2",
    "OG1": " OG1",
    "OG": " OG ",
    "OD2": " OD2",
    "OD1": " OD1",
    "SG": " SG ",
    "SD": " SD ",
    "CG1": " CG1",
    "CG2": " CG2",
    "CG": " CG ",
    "CB": " CB ",

    "CD": " CD ",
    "CD1": " CD1",
    "CD2": " CD2",
    "N": " N  ",
    "C": " C  ",
    "O": " O  ",
    "backbone": frozenset([ " N  ", " CA ", " C  ", " O  " ]),
    "mainchain": frozenset([ " N  ", " CA ", " C  " ])
}


def numsele( string ):
    """ Valid string examples
        :A
        32
        :A.CA
        32:A.CA
        32-40:A
        32.CA
    """
    if isinstance( string, dict ):
        return string
    sele = { "chain": None, "resno": None, "atomname": None }
    if string in [ "*", "", "all" ]:
        return sele
    atomname = string.split(".")
    if len(atomname) > 1 and atomname[1]:
        if len(atomname[1]) > 4:
            raise Exception("atomname must be one to four characters")
        sele["atomname"] = atomname[1][0:4]
    chain = atomname[0].split(":")
    if len(chain) > 1 and chain[1]:
        if len(chain[1]) > 1:
            raise Exception("chain identifier must be one character")
        sele["chain"] = chain[1][0]
    if chain[0]:
        resno = map( int, chain[0].split("-") )
        sele["resno"] = resno[0] if len(resno) == 1 else resno
    return sele


def numdist( numa1, numa2 ):
    return mag( numa1.center() - numa2.center() )


def superpose( npdb1, npdb2, sele1, sele2, subset="CA", inplace=True,
               rmsd_cutoff=None, max_cycles=None, align=True, verbose=True ):

    sele1 = listify( sele1 )
    sele2 = listify( sele2 )

    _sele1 = np.zeros( npdb1.length, bool )
    _sele2 = np.zeros( npdb2.length, bool )

    for s in sele1:
        _sele1 |= npdb1.sele( **s )
    for s in sele2:
        _sele2 |= npdb2.sele( **s )

    numa1 = npdb1.copy( **{ "sele": _sele1, "atomname": subset } )
    numa2 = npdb2.copy( **{ "sele": _sele2, "atomname": subset } )

    msg = []

    def p( m, echo=True ):
        msg.append( m )
        if echo:
            print m

    if align:
        ali1, ali2 = aligner(
            numa1.sequence(), numa2.sequence(),
            method="global", matrix="BLOSUM62"
        )
        ali_sele1 = np.ones( numa1.length, bool )
        ali_sele2 = np.ones( numa2.length, bool )
        i = 0
        j = 0
        for x, y in zip( ali1, ali2 ):
            _i = 0
            _j = 0
            if x == "-":
                ali_sele2[j] = False
            else:
                _i = 1
            if y == "-":
                ali_sele1[i] = False
            else:
                _j = 1
            i += _i
            j += _j
        numa1 = numa1.copy( sele=ali_sele1 )
        numa2 = numa2.copy( sele=ali_sele2 )

    coords1 = numa1['xyz']
    coords2 = numa2['xyz']

    if len( coords1 ) != len( coords2 ):
        raise Exception( "length differ, cannot superpose" )

    if rmsd_cutoff and max_cycles:
        for cycle in xrange( max_cycles ):
            sp = Superposition( coords1, coords2 )
            if verbose:
                p( "Cycle %i, #%i, RMSD %6.2f" % (cycle, sp.n, sp.rmsd) )
            if sp.rmsd <= rmsd_cutoff or sp.n <= 8:
                break
            coords1_trans = sp.transform( coords1, inplace=False )
            d = coords1_trans - coords2
            deviation = np.sqrt( np.sum( d * d, axis=1 ) )
            deviation_idx = np.argsort( deviation )
            coords1 = coords1[ deviation_idx[:-5] ]
            coords2 = coords2[ deviation_idx[:-5] ]
        if cycle == max_cycles - 1:
            p( "warning: still larger than rmsd_cutoff" )
    else:
        sp = Superposition( coords1, coords2 )
    pos = sp.transform( npdb1['xyz'], inplace=inplace )
    npdb1['xyz'] = pos
    numa1['xyz'] = sp.transform( numa1['xyz'], inplace=inplace )
    if verbose:
        p( "RMSD: %f" % sp.rmsd )

    return sp, msg


PDB_DELIMITER = (6, 5, 1, 4, 1, 3, 1, 1, 4, 1, 3, 8, 8, 8, 6, 6, 6, 4, 2, 2 )
PDB_DTYPE = [
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
PDB_USECOLS = (0, 1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 18 )
PDB_COLS = (
    #  0        1         2         3         4         5         6
    ( 0, 6), ( 6, 11), (11, 12), (12, 16), (16, 17), (17, 20), (20, 21),
    #  7         8         9         10        11        12        13
    (21, 22), (22, 26), (26, 27), (27, 30), (30, 38), (38, 46), (46, 54),
    #  14        15        16        17        18        19
    (54, 60), (60, 66), (66, 72), (72, 76), (76, 78), (78, 80)
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
PDB_ATOM_TMPL = (
    "{record:<6}{atomno:>5} {atomname:>4}{altloc:>1}{resname:>3} "
    "{chain:>1}{resno:>4}{insertion:>1}   "
    "{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{bfac:6.2f}      "
    "{segment:<4}{element:>2}{charge:>2}\n"
)

PDB_DEFAULTS2 = {
    "record": "ATOM", "atomno": 1, "atomname": " C  ", "altloc": "",
    "resname": "GLY", "chain": "", "resno": 1, "insertion": "",
    "x": 0.0, "y": 0.0, "z": 0.0,
}
PDB_ATOM_TMPL2 = (
    "{record:<6}{atomno:>5} {atomname:>4}{altloc:>1}{resname:>3} "
    "{chain:>1}{resno:>4}{insertion:>1}   {x:8.3f}{y:8.3f}{z:8.3f}"
    "                          \n"
)

PDB_DEFAULTS3 = {
    "atomno": 1, "atomname": " C  ",
    "resname": "GLY", "chain": "", "resno": 1,
    "x": 0.0, "y": 0.0, "z": 0.0,
}
PDB_ATOM_TMPL3 = (
    "ATOM  {atomno:>5} {atomname:>4} {resname:>3} {chain:>1}{resno:>4}    "
    "{x:8.3f}{y:8.3f}{z:8.3f}                          \n"
)


def numdefaults( natom, defaults ):
    d = {}
    for k, v in defaults.iteritems():
        try:
            d[k] = natom[k]
        except:
            d[k] = v
    return d


def pdb_line( natom, tpl=PDB_ATOM_TMPL, defaults=PDB_DEFAULTS ):
    return tpl.format( **numdefaults( natom, defaults ) )


def xyzr_line( natom ):
    pass


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


SstrucRecord = collections.namedtuple( 'SstrucRecord', [
    'type', 'subtype',
    'chain1', 'resno1', 'resname1',
    'chain2', 'resno2', 'resname2',
    'hbond'
])


class SstrucParser( SimpleParser ):
    def __init__( self ):
        self._list = []

    def _test_line( self, line ):
        return line.startswith("HELIX") or line.startswith("SHEET")

    def _parse_line( self, line ):
        # ensure that res1 comes before res2
        if line.startswith("HELIX"):
            # chain, resno, resname
            res1 = ( line[19], try_int( line[21:25] ), line[15:18] )
            res2 = ( line[31], try_int( line[33:37] ), line[27:30] )
            if res1 > res2:
                res1, res2 = res2, res1
            return SstrucRecord(
                HELIX,
                try_int( line[38:40] ),     # subtype
                res1[0], res1[1], res1[2],  # res1
                res2[0], res2[1], res2[2],  # res2
                None                        # padding...
            )
        elif line.startswith("SHEET"):
            # chain, resno, resname
            res1 = ( line[21], try_int( line[22:26] ), line[17:20] )
            res2 = ( line[32], try_int( line[33:37] ), line[28:31] )
            if res1 > res2:
                res1, res2 = res2, res1
            return SstrucRecord(
                SHEET,
                try_int( line[38:40] ),     # strand sense (subtype)
                res1[0], res1[1], res1[2],  # res1
                res2[0], res2[1], res2[2],  # res2
                try_int( line[65:69], False ),     # resno hbond prev strand
            )

    def get( self ):
        self._list.sort( key=operator.attrgetter("chain1", "resno1") )
        # for ss in self._list:
        #     print ss
        filtered = []
        it = iter_window( self._list, 2 )
        for ss1, ss2 in it:
            # remove sstruc that occur twice i.e. in 3s0c
            # once as the starting strand and then making an
            # hbond - take the second one
            if (ss1.chain1, ss1.resno1) == (ss2.chain1, ss2.resno1):
                continue
            # remove sstruc that is within another i.e. in 2plv
            # take the one that starts first
            # NO overlaps are OK and handled later
            # if (ss1.chain2, ss1.resno2)>(ss2.chain1, ss2.resno1):
            #     next(it)
            filtered.append( ss1 )
        self._list = filtered
        return self._list

MissingRecord = collections.namedtuple( 'MissingRecord', [
    'type', 'model',
    'resname', 'chain', 'ssseqi',
    'identifier', 'atoms'
])


class MissingParser( SimpleParser ):
    def __init__( self ):
        self._list = []

    def _test_line( self, line ):
        return line.startswith("REMARK 465") or line.startswith("REMARK 470")

    def _parse_line( self, line ):
        if line.startswith("REMARK 465"):
            model = line[11:15]
            resname = line[15:18]
            chain = line[18:20]
            ssseq = line[20:26]
            identifier = line[26:28]
            if try_int(ssseq, 'False' ) != 'False':
                return MissingRecord(
                    "Residue",
                    model,     # model
                    resname,  # resname
                    chain,  # chain
                    try_int(ssseq, False ),  # ssseqi
                    identifier,
                    None  # atoms
                )
            #else:
                #return MissingRecord(None,None,None,None,None,None,None)
        elif line.startswith("REMARK 470"):
            model = line[13:15]
            resname = line[15:18]
            chain = line[19]
            ssseq = line[20:25]
            identifier = line[25:28]
            atoms = line[28:].rstrip()
            if try_int(ssseq, 'False' ) != 'False':
                return MissingRecord(
                    "Atoms",
                    model,     # model
                    resname,  # resname
                    chain,  # chain
                    try_int(ssseq, False ),  # ssseqi
                    identifier,
                    atoms  # atoms
                )
            #else:
                #return MissingRecord(None,None,None,None,None,None,None)

    def get( self ):
        return self._list


class InfoParser( object ):
    def __init__( self ):
        self._dict = collections.defaultdict( str )

    def __call__( self, line ):
        self._parse_line( line )

    def _parse_line( self, line ):

        if line.startswith("KEYWDS"):
            self._dict["keywords"] += line[10:].rstrip()
        elif line.startswith("EXPDTA"):
            self._dict["experiment"] += line[10:].rstrip()
        elif line.startswith("MDLTYP"):
            self._dict["model_type"] += line[10:].rstrip()
        elif line.startswith("TITLE"):
            self._dict["title"] += line[10:].rstrip()
        elif line.startswith("SPLIT"):
            self._dict["split"] += line[10:].rstrip()
        elif line.startswith("REMARK"):
            if line[9] == "2":
                self._dict["resolution"] += line[10:].rstrip()
        elif line.startswith("HEADER"):
            if self._dict["header"]:
                self._dict["header"].append(line[10:50].rstrip()+ line[50:62].rstrip()+ line[62:].rstrip())
            else:
                self._dict["header"] = [line[10:50].rstrip(), line[50:62].rstrip(), line[62:].rstrip()]
        elif line.startswith("OBSLTE"):
            self._dict["obsolete"] += line[31:].rstrip()

    def get( self ):
        dct = self._dict
        mdl = [ s.strip() for s in dct.get("model_type", "").split(";") ]
        mdl_dct = {}
        for typ in mdl:
            x = s.split(",", 1)
            if len(x) == 2:
                mdl_dct[ x[0].strip() ] = x[1].strip()
            elif len(x) == 1 and x[0]:
                mdl_dct[ x[0].strip() ] = True
        return {
            "keywords": dct.get("keywords", "").replace("- ", "-").split(", "),
            "experiment": dct.get("experiment", ""),
            "title": dct.get("title", ""),
            "resolution": try_float(
                get_index(
                    dct.get("resolution", "").split(), 1, None
                ),
                None
            ),
            "splited_entry": dct.get("split", "").split(),
            "model_type": mdl_dct,
            "obsolete": dct.get("obsolete", "").split(),
            "header": dct.get("header", "")
        }


# TODO needs documentation
_BORDER = [{
    'resno': None, 'insertion': None, 'chain': None,
    'sstruc': None, 'altloc': None
}]


def BORDER( _atoms ):
    return itertools.chain( _atoms, _BORDER )


class NumAtoms:
    def __init__( self, atoms, coords, flag=None ):
        self._atoms = atoms
        self._coords = coords
        self.flag = flag
        self.length = len( atoms )
        self._memo = {}

    def __len__( self ):
        return len( self._atoms )

    def __getitem__( self, key ):
        if key == 'xyz':
            return self._coords
        else:
            return self._atoms[ key ]

    def __setitem__( self, key, value ):
        if key == 'xyz':
            for i, d in enumerate('xyz'):
                self._atoms[ d ] = value[ ..., i ]
            self._coords = np.vstack((
                self._atoms['x'], self._atoms['y'], self._atoms['z']
            )).T
            # print "may_share_memory: %s" % np.may_share_memory( self._atoms, self._coords )
            # print "flags.owndata: %s, %s" % (
            #     self._atoms.flags.owndata, self._coords.flags.owndata
            # )
        else:
            self._atoms[ key ] = value

    def sele( self, chain=None, resno=None, atomname=None, atomno=None,
              altloc=None, resname=None, sele=None, record=None,
              invert=None ):
        atoms = self._atoms
        if sele is None:
            sele = np.ones( self.length, bool )
        if record is not None:
            sele &= atoms['record'] == record
        if chain is not None:
            if isinstance( chain, collections.Iterable ):
                tmps = atoms['chain'] == chain[0]
                for an in chain[1:]:
                    tmps |= atoms['chain'] == an
                sele &= tmps
            else:
                sele &= atoms['chain'] == chain
        if resname is not None:
            if ( isinstance( resname, collections.Iterable )
                    and not isinstance( resname, basestring ) ):
                print resname
                tmps = atoms['resname'] == resname[0]
                for an in resname[1:]:
                    tmps |= atoms['resname'] == an
                sele &= tmps
            else:
                sele &= atoms['resname'] == resname
        if resno is not None:
            if isinstance( resno, collections.Iterable ):
                sele &= (
                    (atoms['resno'] >= resno[0]) &
                    (atoms['resno'] <= resno[1])
                )
            else:
                sele &= atoms['resno'] == resno
        if atomno is not None:
            if isinstance( atomno, collections.Iterable ):
                sele &= (
                    (atoms['atomno'] >= atomno[0]) &
                    (atoms['atomno'] <= atomno[1])
                )
            else:
                sele &= atoms['atomno'] == atomno
        if atomname is not None:
            if ( isinstance( atomname, collections.Iterable )
                    and not isinstance( atomname, basestring ) ):
                atomname = [ ATOMS.get( a, a ) for a in atomname ]
                tmps = atoms['atomname'] == atomname[0]
                for an in atomname[1:]:
                    tmps |= atoms['atomname'] == an
                sele &= tmps
            else:
                atomname = ATOMS.get( atomname, atomname )
                sele &= atoms['atomname'] == atomname
        if altloc is not None:
            if isinstance( altloc, collections.Iterable ):
                tmps = atoms['altloc'] == altloc[0]
                for an in altloc[1:]:
                    tmps |= atoms['altloc'] == an
                sele &= tmps
            else:
                sele &= atoms['altloc'] == altloc
        if invert:
            np.logical_not( sele, sele )
        return sele

    def slice( self, begin, end=None, flag=None ):
        if not end:
            end = begin + 1
        return NumAtoms(
            self._atoms[begin:end],
            self._coords[begin:end], flag=flag
        )

    def copy( self, **sele ):
        coords, atoms = self._select( **sele )
        return NumAtoms( atoms, coords )

    def _select( self, **sele ):
        coords = self._coords
        atoms = self._atoms
        if len(sele):
            _sele = self.sele( **sele )
            if len(_sele):
                coords = self._coords[ _sele ]
                atoms = self._atoms[ _sele ]
            else:
                coords = np.array( [], dtype=coords.dtype )
                atoms = np.array( [], dtype=atoms.dtype )
        return coords, atoms

    def get( self, key, **sele ):
        coords, atoms = self._select( **sele )
        if key == 'xyz':
            return coords
        else:
            return atoms[ key ]

    def index( self, first=False, last=False, interval=False, **sele ):
        indices = np.nonzero( self.sele( **sele ) )[0]
        if first:
            return indices[0]
        elif last:
            return indices[-1]
        elif interval:
            return indices[0], indices[-1]
        else:
            return indices

    @memoize_m
    def iter_chain( self, **sele ):
        memo = []
        coords, atoms = self._select( **sele )
        if len(atoms) > 0:
            chain = atoms['chain'][0]
            k = 0
            l = 0
            for a in BORDER( atoms ):
                if chain != a['chain']:
                    numa = self.slice( k, l )
                    memo.append( numa )
                    chain = a['chain']
                    k = l
                l += 1
        return memo

    @memoize_m
    def iter_sstruc( self, **sele ):
        memo = []
        for numatoms in self.iter_chain( **sele ):
            sstruc = numatoms['sstruc'][0]
            atoms = numatoms._atoms
            k = 0
            l = 0
            for a in BORDER( atoms ):
                if sstruc != a['sstruc']:
                    memo.append( numatoms.slice( k, l ) )
                    k = l
                    sstruc = a['sstruc']
                l += 1
        return memo

    @memoize_m
    def _iter_resno( self, **sele ):
        memo = []
        for numatoms in self.iter_chain( **sele ):
            res = (
                numatoms['resno'][0],
                numatoms['insertion'][0],
                numatoms['altloc'][0]
            )
            k = 0
            l = 0
            for a in BORDER( numatoms._atoms ):
                if res != ( a['resno'], a['insertion'], a['altloc'] ):
                    memo.append( numatoms.slice( k, l ) )
                    k = l
                    res = ( a['resno'], a['insertion'], a['altloc'] )
                l += 1
        return memo

    @memoize_m
    def iter_resno( self, incomplete=False, **sele ):
        memo = []
        numa_prev = None
        numa0_prev = None
        for numa in self._iter_resno( **sele ):
            numa0 = numa[0]
            if incomplete or not numa0['incomplete']:
                numa.flag = True
                if numa_prev:
                    # skip any additional altloc
                    if numa0['resno'] == numa0_prev['resno'] and \
                            numa0['insertion'] == numa0_prev['insertion'] and \
                            numa0['altloc'] != numa0_prev['altloc']:
                        continue
                    if numa0_prev['incomplete'] or numa0['incomplete']:
                        numa.flag = True
                    elif numa0_prev['resno'] == numa0['resno'] - 1:
                        numa.flag = False
                    elif numdist(
                            # assuming the first three atoms are N, CA, C
                            numa_prev.slice( 2 ),
                            numa.slice( 0 ) ) > 1.4:
                            # numa_prev.copy( atomname="C" ),
                            # numa.copy( atomname="N" ) ) > 1.4:
                        numa.flag = True
                    else:
                        numa.flag = False
                memo.append( numa )
                numa_prev = numa
                numa0_prev = numa0
        return memo

    @memoize_m
    def iter_resno2( self, window, **sele ):
        memo = []
        # (TODO) assumes the first a has a.flag==True
        it = iter( self.iter_resno( **sele ) )
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
                    break
            else:
                result = result[1:] + (a,)
            memo.append( result )
        return memo

    def axis( self, **sele ):
        try:
            return axis( self.get( 'xyz', **sele ) )
        except Exception as e:
            LOG.error( "axis (%s) => %s, %s" % ( e, sele, self._atoms ) )
            raise e

    def sequence( self, as_array=False, **sele ):
        s = [
            AA1.get( a['resname'][0], "?" )
            for a in self._iter_resno( **sele )
        ]
        if as_array:
            return np.array( s, dtype="|S1" )
        else:
            return "".join( s )

    def center( self, **sele ):
        coords = self.get( 'xyz', **sele )
        try:
            return np.sum( coords, axis=0 ) / len(coords)
        except FloatingPointError as e:
            print self._atoms
            raise e

    def dist( self, sele1, sele2 ):
        return mag( self.center( **sele1 ) - self.center( **sele2 ) )

    def write( self, file_name, order='original',
               tpl=PDB_ATOM_TMPL, defaults=PDB_DEFAULTS, **sele ):
        coords, atoms = self._select( **sele )
        if order == 'original':
            atoms = np.sort( atoms, order='atomno' )
        with open( file_name, "wb" ) as fp:
            for natom in atoms:
                fp.write( pdb_line( natom, tpl=tpl, defaults=defaults ) )
            fp.write( "END" )
            # for numa in self.iter_resno( **sele ):
            #     for a in numa._atoms:
            #         fp.write( pdb_line( a ) )

    def write2( self, file_name, **sele ):
        """ write only the first conformation and
            set altloc empty
        """
        coords, atoms = self._select( **sele )
        na = NumAtoms( atoms, coords )
        with open( file_name, "wb" ) as fp:
            i = 1
            for numa in na.iter_resno( incomplete=True ):
                for a in numa._atoms:
                    a["altloc"] = " "
                    a["atomno"] = i
                    fp.write( pdb_line( a ) )
                    i += 1
            fp.write( "END" )


class NumPdb:
    numatoms = None

    def __init__( self, pdb_path, features=None ):
        self.pdb_path = pdb_path
        self.pdb_id = basekit.utils.path.stem( pdb_path )
        self.features = {
            "phi_psi": True,
            "sstruc": True,
            "detect_missing": False,
            "backbone_only": False,
            "protein_only": True,
            "detect_incomplete": True,
            "configuration": True,
            "no_sort": False,
            "info": False,
            "extra": False
        }
        if features:
            self.features.update( features )
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
        if self.features["configuration"]:
            types += [ ( 'configuration', np.bool ) ]
            extra += [ False ]
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
        if self.features["info"]:
            parsers[ "info" ] = InfoParser()
        if self.features["detect_missing"]:
            parsers[ "missing" ] = MissingParser()
        if self.features["extra"]:
            for xname, xtype, xdefault in self.features["extra"]:
                types += [ ( xname, xtype ) ]
                extra += [ xdefault ]

        atoms = []
        atoms_append = atoms.append
        header = []
        header_append = header.append

        with open( self.pdb_path, "rb" ) as fp:
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
                if line[0:4] == "ATOM" or \
                        ( line[0:6] == "HETATM" and line[17:20] in aminoacids ):
                    if po and line[17:20] in nucleotides:
                        continue
                    if ( nbo or line[12:16] in backbone ) and line[16] in altloc:
                        atoms_append( tupl(
                            [ line[ c[0]:c[1] ] for c in cols ] + extra
                        ))
                elif line[0:6] == "HETATM":
                    atoms_append( tupl(
                        [ line[ c[0]:c[1] ] for c in cols ] + extra
                    ))
                elif line[0:5] == "MODEL":
                    # TODO add model field
                    pass
                elif line[0:6] == "ENDMDL":
                    # stop condition - for now
                    break
                elif line[0:6] == "CONECT":
                    # stop condition
                    break

        for name in parsers.keys():
            p = parsers.get( name )
            if p:
                self.__dict__[ "_%s" % name ] = p.get()

        atoms = np.array(atoms, dtype=types)
        if not self.features["no_sort"]:
            atoms.sort( order=[
                'chain', 'resno', 'insertion', 'altloc'
            ])
        coords = np.vstack(( atoms['x'], atoms['y'], atoms['z'] )).T
        # print "may_share_memory: %s" % np.may_share_memory( atoms, coords )
        # print "flags.owndata: %s, %s" % (
        #     atoms.flags.owndata, coords.flags.owndata
        # )

        self._atoms = atoms
        self._coords = coords
        self._header = header
        self.numatoms = NumAtoms( self._atoms, self._coords )
        self.length = len( atoms )

        if self.features["configuration"]:
            self.__calc_configuration()
        if self.features["detect_incomplete"]:
            self.__calc_incomplete()
        if self.features["phi_psi"]:
            self.__calc_phi_psi()
        if self.features["sstruc"]:
            self.__calc_sstruc()

    def __calc_configuration( self ):
        pass

    def __calc_incomplete( self ):
        bb_subset = ATOMS['backbone'].issubset
        try:
            numa = False
            for numa in self._iter_resno():
                if not bb_subset( numa["atomname"] ):
                    numa['incomplete'] = True
                    #print numa._atoms
        except Exception as e:
            LOG.error( "[%s] calc incomplete (%s) => %s" % (
                self.pdb_id, e, numa and numa['atomno'][0]
            ))

    def __calc_sstruc( self ):
        for ss in self._sstruc:
            try:
                idx_beg = self.index(
                    chain=ss.chain1, resno=ss.resno1, first=True
                )
                idx_end = self.index(
                    chain=ss.chain2, resno=ss.resno2, last=True
                )
                if ss.type == HELIX:
                    self.slice( idx_beg, idx_end + 1 )['sstruc'] = "H"
                elif ss.type == SHEET:
                    self.slice( idx_beg, idx_end + 1 )['sstruc'] = "E"
            except Exception as e:
                LOG.error( "[%s] calc sstruc (%s) => %s" % (
                    self.pdb_id, e, ss
                ))

    def __calc_phi_psi( self, dihedral=dihedral ):
        mainchain = ATOMS['mainchain']
        for na_curr, na_next in self.iter_resno2( 2 ):
            try:
                # xyz_n, xyz_ca, xyz_c = na_curr.get(
                #     'xyz', atomname=mainchain
                # )
                # xyz_n_next, xyz_ca_next, xyz_c_next = na_next.get(
                #     'xyz', atomname=mainchain
                # )
                # assuming the first three atoms are N, CA, C
                xyz_n, xyz_ca, xyz_c = na_curr._coords[0:3]
                xyz_n_next, xyz_ca_next, xyz_c_next = na_next._coords[0:3]
                na_curr['psi'] = dihedral(
                    xyz_n, xyz_ca, xyz_c, xyz_n_next
                )
                na_next['phi'] = dihedral(
                    xyz_c, xyz_n_next, xyz_ca_next, xyz_c_next
                )
            except Exception as e:
                LOG.error( "[%s] calc phi/psi (%s) => %s, %s, %s, %s, %s, %s" % (
                    self.pdb_id, e,
                    na_curr['resno'], na_curr['chain'],
                    na_next['resno'], na_next['chain'],
                    na_curr['atomname'], na_next['atomname']
                ))
        # for a in self.iter_resno(): print a._atoms[0], a.flag
