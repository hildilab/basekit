#! /usr/bin/env python
from __future__ import division

import argparse
import os
import os.path
from collections import defaultdict
import logging
import re
from string import Template

HOLE_PARTLY_FILLED = 1
HOLE_NOT_FILLED = 2
HOLE_PARTLY_FILLED_HETS_REMOVED = 3
HOLE_FILLED_HETS_REMOVED = 4

logging.basicConfig()
LOG = logging.getLogger('prep')
LOG.setLevel( logging.ERROR )


def get_index(seq, index, default=None):
    if not seq:
        return default
    try:
        return seq[index]
    except IndexError:
        return default


def memoize(f):
    cache = {}
    def memf(*x):
        if x not in cache:
            cache[x] = f(*x)
        return cache[x]
    return memf


@memoize
def get_pdb_coord_dict(pdb_file):
    pdb_fp = open(pdb_file)
    coord_dict = {}
    i = 1
    for l in pdb_fp:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            key = ( float(l[30:38]), float(l[38:46]), float(l[46:54]) )
            if key in coord_dict:
                LOG.error( "coords already in dict. %s" % str(key) )
            else:
                coord_dict[ key ] = i
                i += 1
    pdb_fp.close()
    return coord_dict


@memoize
def get_pdb_index_dict(pdb_file):
    pdb_fp = open(pdb_file)
    index_dict = {}
    i = 1
    for l in pdb_fp:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            index_dict[ i ] = l
            i += 1
    pdb_fp.close()
    return index_dict


def prep_contact(contact_file, pdb_file):
    contact_fp = open(contact_file)
    name, ext = os.path.splitext(contact_file)
    elms_out_fp = open( "%s.atmsele" % (contact_file), "w")
    contact_out_fp = open( "%s.atmprop" % (contact_file), "w")
    contact2_out_fp = open( "%s_contact.atmsele" % (contact_file), "w")

    pdb_coord_dict = get_pdb_coord_dict(pdb_file)
    pdb_index_dict = get_pdb_index_dict(pdb_file)
    contact_index_dict = {}
    structure_elms = []
    contact_lines = [None] * len(pdb_coord_dict)

    i = 1
    for l in contact_fp:
        if l.startswith('OVERVW') or l.startswith('OVERV2'):
            structure_elms += l.split()[1:]
        elif l.split(' ')[0] in ['END', 'CUTOFF', 'HEADER', 'COMPND', 'REMARK', 'METHOD', 'LENGTH']:
            pass
        elif l[0] in ['H', 'C', 'E', 'O', 'W']:
            key = ( float(l[27:35]), float(l[35:43]), float(l[43:51]) )
            if key in pdb_coord_dict:
                j = pdb_coord_dict[key]
                contact_lines[ j-1 ] = l
                contact_index_dict[ i ] = j
                i += 1
            else:
                LOG.error( "contact coords not in pdb coords dict. %s" % str(key) )

    if len(pdb_coord_dict) != len(contact_index_dict):
        LOG.warning( "number of atoms in contact and pdb coord dicts differs." )

    structure_elms_index = {}
    for i, elm in enumerate(structure_elms):
        structure_elms_index[ i+1 ] = elm
    structure_elms_index[ 0 ] = "Water"
    structure_elms_index[ -1 ] = "Membrane"
    structure_elms_index[ -51 ] = "Membrane"
    # -50, 50, 67

    elms_list = ['Membrane', 'Water'] + structure_elms
    elms_list_dict = {}
    for i, e in enumerate(elms_list):
        elms_list_dict[e] = i

    structure_elms_dict = defaultdict(list)
    contacts_cutoff_dict = defaultdict(list)

    contact_out_fp.write('%s\n' % " ".join([ "%s#9" % e for e in elms_list ]))
    for i, l in enumerate(contact_lines):
        elms_out_list = ['9'] * len(elms_list)
        if l:
            contacts = l[73:].split()
            if len(contacts) % 2 != 0:
                LOG.error( "contacts list has uneven number of elements. %s" % l )
            for elm_id, cutoff in zip(contacts[::2],contacts[1::2]):
                elm_id = int(elm_id)
                if elm_id in structure_elms_index:
                    elm_name = structure_elms_index[ int(elm_id) ]
                    elms_out_list[ elms_list_dict[ elm_name ] ] = cutoff
                    # zero based atomindex
                    contacts_cutoff_dict[ "%s_%s" % (elm_name, cutoff) ].append( i )
                else:
                    # LOG.warning( "structure element not known. %s" % elm_id )
                    LOG.debug( "structure element not known. %s" % elm_id )
            structure_elms_dict[ l[0:5] ].append( i+1 )
        else:
            LOG.warning( "no contact data for line (@%s): %s" % (i, pdb_index_dict[ i+1 ].strip('\n')) )
        contact_out_fp.write( '%s\n' % " ".join( elms_out_list ) )

    for name, atoms in structure_elms_dict.iteritems():
        # zero based atomindex
        struc = " ".join([ str(a-1) for a in atoms ])
        elms_out_fp.write( "%s %s\n" % (name, struc) )

    for name, atoms in contacts_cutoff_dict.iteritems():
        contacts = " ".join([ str(a) for a in atoms ])
        contact2_out_fp.write( "%s %s\n" % (name, contacts) )
        
    contact_fp.close()
    elms_out_fp.close()
    contact_out_fp.close()
    contact2_out_fp.close()




def prep_volume(vol_file, pdb_file):
    vol_fp = open(vol_file)
    name, ext = os.path.splitext(vol_file)
    holes_out_fp = open( "%s.atmsele" % (vol_file), "w")
    prop_out_fp = open( "%s.atmprop" % (vol_file), "w")

    pdb_coord_dict = get_pdb_coord_dict(pdb_file)
    pdb_index_dict = get_pdb_index_dict(pdb_file)
    vol_index_dict = {}

    vol_lines = [None] * len(pdb_coord_dict)
    hole_lines = []
    nrholes = {}
    hole_types = []
    
    i = 1
    for l in vol_fp:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            key = ( float(l[30:38]), float(l[38:46]), float(l[46:54]) )
            if key in pdb_coord_dict:
                j = pdb_coord_dict[key]
                vol_lines[ j-1 ] = l
                vol_index_dict[ i ] = j
                i += 1
            else:
                LOG.error( "vol coords not in pdb coords dict. %s" % str(key) )
        elif l.startswith('HOLE NUMBER'):
            hole_lines.append( l )
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
        LOG.warning( "number of atoms in vol and pdb coord dicts differs." )

    for l in hole_lines:
        ls = l[12:].split()
        neighbours = []
        for nb in ls[1:]:
            nb = int(nb)
            if nb in vol_index_dict:
                neighbours.append( vol_index_dict[nb] )
            else:
                LOG.error( "hole neighbour index not found. %s" % nb )
        neighbours.sort()
        # zero based atomindex
        neighbours = " ".join([ str(nb-1) for nb in neighbours ])
        holes_out_fp.write( "HOLE_NUMBER_%s_%i %s\n" % (
            ls[0], hole_types[ int(ls[0])-1 ], neighbours
        ))

    prop_out_fp.write('volume_vdw#-1 volume_vdw_1_4#-1 buried_flag#-1 packing_density#-1\n')
    for i, l in enumerate(vol_lines):
        if l:
            ls = l.split()
            vdwvol = float(ls[-3])  #VOLUME INSIDE VAN-DER-WAALS SPHERE
            sevol = float(ls[-2])   #VOLUME IN 1.4 ANGSTROM LAYER OUTSIDE VDW-SPHERE
            if vdwvol==0.0:
                packdens = 0.0
                LOG.error( "vdw volume zero. %s" % l )
            elif (vdwvol+sevol)==0.0:
                packdens = 999.99
                LOG.error( "sum of vdw volume and excluded volume zero. %s" % l )
            else:
                packdens = (vdwvol/(vdwvol+sevol))
            out_line = '%s %.2f\n' % (l[67:82], packdens)
        else:
            out_line = '-1 -1 -1 -1\n'
            LOG.warning( "no volume data for line: %s" % pdb_index_dict[ i+1 ].strip('\n') )
        prop_out_fp.write( out_line )

    vol_fp.close()
    prop_out_fp.close()
    holes_out_fp.close()



def prep_pdb(pdb_file):
    pdb_fp = open(pdb_file)
    name, ext = os.path.splitext(pdb_file)
    pdb_out_fp = open( "%s_short%s" % (name, ext), "w")
    
    for line in pdb_fp:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            pdb_out_fp.write( line[:54]+'\n' )

    pdb_fp.close()
    pdb_out_fp.close()


@memoize
def get_pdb_chainres_dict(pdb_file):
    pdb_fp = open(pdb_file)
    chainres_dict = defaultdict(list)
    i = 1
    for l in pdb_fp:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            key = ( l[21], int(l[22:26]) )
            chainres_dict[ key ].append( i )
            i += 1
    pdb_fp.close()
    return chainres_dict



def prep_tmhelix(tmhelix_file, pdb_file):
    tmhelix_fp = open(tmhelix_file)
    tmhelix_out_fp = open( "%s.atmsele" % (tmhelix_file), "w")

    pdb_chainres_dict = get_pdb_chainres_dict(pdb_file)

    i = 1
    for l in tmhelix_fp:
        ls = l.split()
        if len(ls) == 3:
            # the chain identifier can be blank ' '
            keys = [ ( ' ', int(ls[1]) ), ( ' ', int(ls[2]) ) ]
        elif len(ls) == 5:
            keys = [ ( ls[1], int(ls[2]) ), ( ls[3], int(ls[4]) ) ]
        else:
            continue
        indexes = []
        for key in keys:
            if key in pdb_chainres_dict:
                indexes.append( pdb_chainres_dict[ key ] )
            else:
                LOG.error( "key not found in chainres dict. %s" % key )
        # zero based atomindex
        tmhelix_out_fp.write( "TMHELIX_%s %s:%s\n" % (i, indexes[0][0]-1, indexes[1][-1]-1) )
        i += 1
    tmhelix_fp.close()
    tmhelix_out_fp.close()


@memoize
def get_pdb_chainresatom_dict(pdb_file):
    pdb_fp = open(pdb_file)
    chainresatom_dict = {}
    i = 1
    for l in pdb_fp:
        if l.startswith('ATOM') or l.startswith('HETATM'):
            key = ( l[21], int(l[22:26]), l[12:16].strip() )
            if key in chainresatom_dict:
                LOG.error( "key already in chainresatom dict. %s" % str(key) )
            else:
                chainresatom_dict[ key ] = i
            i += 1
    pdb_fp.close()
    return chainresatom_dict


def prep_hbexplore(hbexplore_file, pdb_file):
    hbexplore_fp = open(hbexplore_file)
    hbexplore_out_fp = open( "%s.bonds" % (hbexplore_file), "w")
    hbexplore_out_fp.write('data "connect_atoms"\n')

    chainresatom_dict = get_pdb_chainresatom_dict(pdb_file)
    id_dict = {}
    hbpart=0
    #   SG  CYS   127      O   ILE   124    3.2597     
    # TOM     18  CB  VAL     2       3.462  14.126   6.256  1.00 14.07      193L 162
    # parse HBX.anal file
    #87     NH2 ARG A  71      O   MET A   1    3.0821  2.356   128.0    53.8   s-s
    for l in hbexplore_fp:
        if not l:
            continue
        if l.startswith('No.'):
            hbpart=1
        elif hbpart and l[0]=='=':
            break
        elif not hbpart or len(l)<=30:
            continue
        elif hbpart and l[0]!='-':
            # print l[:-1]
            ts = l.split("\t")
            t=ts[1]
            id = int(ts[0])
            if( id not in id_dict ):
                id_dict[ id ] = True
                # key = ( chain, resno, atomname )
                key1 = (
                    t[12:14].strip(),   # chain
                    int(t[14:18]),      # resno
                    t[4:7].strip()      # atomname
                )
                key2 = (
                    t[31:33].strip(),   # chain
                    int(t[33:37]),      # resno
                    t[23:27].strip()    # atomname
                )
                if key1 in chainresatom_dict:
                    idx1 = chainresatom_dict[ key1 ]
                else:
                    LOG.error( "key not found in chainresatom dict. %s" % str(key1) )
                    continue
                if key2 in chainresatom_dict:
                    idx2 = chainresatom_dict[ key2 ]
                else:
                    LOG.error( "key not found in chainresatom dict. %s" % str(key2) )
                    continue
                # zero based atomindex
                # hbexplore_out_fp.write( '%s %s 2049 0.1 0.0 hbond;\n' % ( idx1-1, idx2-1 ) )
                hbexplore_out_fp.write( '%s %s 2048 0.1 0.0 hbond;\n' % ( idx1-1, idx2-1 ) )
        
    hbexplore_out_fp.write( 'end "connect_atoms";\n' )

    hbexplore_fp.close()
    hbexplore_out_fp.close()


def create_json( values_dict, out, tpl ):
    with open( tpl, "r" ) as fp:
        tpl_str = fp.read()
    with open( out, "w" ) as fp:
        fp.write( Template( tpl_str ).substitute( **values_dict ) )



def main():

    # create the parser
    parser = argparse.ArgumentParser(
        description = __doc__,
    )
    # add the arguments
    parser.add_argument('-o', default='.', help='output directory')
    parser.add_argument('-pdb', help='pdb file')
    parser.add_argument('-vol', help='volume file')
    parser.add_argument('-contact', help='contact (sco, mbn) file')
    parser.add_argument('-tmhelix', help='tmhelix file')
    parser.add_argument('-hbx', help='hbexplore file')
    parser.add_argument('-jsontpl', help='json template file')
    parser.add_argument('-dirwalk', nargs='+', default=[], help='dir, [regex]')

    # parse the command line
    args = parser.parse_args()

    if args.pdb:
        prep_pdb( args.pdb )

    if args.vol and args.pdb:
        prep_volume( args.vol, args.pdb )

    if args.contact and args.pdb:
        prep_contact( args.contact, args.pdb )

    if args.tmhelix and args.pdb:
        prep_tmhelix( args.tmhelix, args.pdb )

    if args.hbx and args.pdb:
        prep_hbexplore( args.hbx, args.pdb )

    if args.dirwalk:
        print args.dirwalk
        directory = args.dirwalk[0]
        pattern = get_index( args.dirwalk, 1, False )
        for pathname in os.listdir( directory ):
            if os.path.isfile( pathname ):
                continue
            match = ''
            if pattern:
                m = re.match( pattern, pathname )
                if m: 
                    match = m.group(1)
                else:
                    continue
            print pathname, match
            if args.jsontpl:
                values_dict = { "id": match }
                out = os.path.join( directory, pathname, "json.provi" )
                create_json( values_dict, out, args.jsontpl )


            


if __name__ == "__main__":
    main()