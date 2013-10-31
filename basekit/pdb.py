from __future__ import with_statement

import os
import re
import urllib2
import gzip
import collections
import json
import string
import datetime

import numpy as np
np.seterr( all="raise" )

import utils.numpdb as numpdb
from utils import try_int, copy_dict
from utils.tool import _, _dir_init, PyTool, ProviMixin
from utils.timer import Timer
from utils.db import get_pdb_files
from utils.listing import ListRecord, ListIO, list_compare, list_join
from utils.math import rmatrixu


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "pdb" )

PDB_SEARCH_URL = 'http://www.rcsb.org/pdb/rest/search'
PDB_DOWNLOAD_URL = 'http://www.rcsb.org/pdb/files/'
PDBE_ASSEMBLY = 'http://www.ebi.ac.uk/pdbe-srv/view/files/{pdb_id:}_1.mmol'

DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "motif" )
BASEKIT_DIR = os.path.split( PARENT_DIR )[0]
DATA_DIR = os.path.join( BASEKIT_DIR, "data", "pdb" )


class SplitPdbSSE (PyTool) :
    args = [
            _( "pdb_file", type="dir" )            
        ]
    out = [
            _( "pieces_dir", dir="pieces" )
        ]
    def func( self ):
        npdb = numpdb.NumPdb( self.pdb_file, {
            "phi_psi": False,
            "sstruc": True,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False,
            "info": False
        })
        for i, numa in enumerate( npdb.iter_sstruc() ):
            print numa["chain"][0], numa["resno"].min(), numa["resno"].max()
            print numa.sequence()
            values_dic ={}
            values_dic[1]=numa['resno'][0]
            values_dic[2]=numa['resno'][-1]
            values_dic[3]=numa.sequence()
            values_dic[4]=numa['chain'][0]
            if not os.path.exists(self.pieces_dir):
                os.makedirs(self.pieces_dir)
            file_name=os.path.join(self.pieces_dir,  "%s_%s_%i-%i.%s" % (self.pdb_file[-8:-4],numa['chain'][0],numa['resno'][0],numa['resno'][-1],'pdb'))
            file_namej=os.path.join(self.pieces_dir, "%s_%s_%i-%i.%s" % (self.pdb_file[-8:-4],numa['chain'][0],numa['resno'][0],numa['resno'][-1],'json'))
            numa.write(file_name)
            with open(file_namej, "w" ) as fp:
                json.dump( values_dic,fp, indent=4 )



def parse_het_dictionary( het_file ):
    het_dict = collections.defaultdict( list )
    key = None
    with open( het_file, "r" ) as fp:
        for line in fp:
            if line.startswith("RESIDUE"):
                key = line[10:13]
            if line.startswith("CONECT"):
                het_dict[ key ].append( line[11:15] )
    print len(het_dict)
    aminoacid_list = []
    for key, atom_list in het_dict.iteritems():
        if numpdb.ATOMS['backbone'].issubset( atom_list ):
            aminoacid_list.append( key )
    print len(aminoacid_list)
    return aminoacid_list
    
class PdbHetDictionary( PyTool ):
    """
    Tool to extract the names of (non-standard) amino acids from
    ftp://ftp.rcsb.org/pub/pdb/data/monomers/het_dictionary.txt
    """
    args = [
        _( "het_file", type="file", ext="txt", help="hetero dictionary" )
    ]
    out = [
        _( "aa_file", file="aa.json" )
    ]
    def func( self ):
        aminoacid_list = parse_het_dictionary( self.het_file )
        with open( self.aa_file, "w" ) as fp:
            json.dump( aminoacid_list, fp )




def unzip_pdb( fpath ):
    fdir, fname = os.path.split( fpath )
    fpdb = os.path.join( fdir, fname[3:7]+".pdb" )
    with open( fpdb, "w" ) as fp:
        with gzip.open( fpath, "rb" ) as zp:
            fp.write( zp.read() )
    return True

class PdbUnzip( PyTool ):
    args = [
        _( "pdb_archive", type="dir" )
    ]
    def func( self ):
        pdb_file_list = get_pdb_files( self.pdb_archive, pattern="ent.gz" )
        for pdb_file in pdb_file_list:
            unzip_pdb( pdb_file )



def pdb_download( pdb_id, output_file ):
    pdb_download_urls = [
        PDB_DOWNLOAD_URL,
        "http://198.202.122.52/pdb/files/", # outdated?
        "http://198.202.122.51/pdb/files/"  # outdated?
    ]
    for base_url in pdb_download_urls:
        try:
            url = "%s%s.pdb" % ( base_url, pdb_id )
            data = urllib2.urlopen( url ).read()
            with open( output_file, 'wb' ) as fp:
                fp.write( data )
            return True
        except urllib2.HTTPError:
            pass
    else:
        raise Exception( "Error downloading pdb entry '%s'" % pdb_id )

class PdbDownload( PyTool ):
    args = [
        _( "pdb_id", type="text", 
            help="single id or multiple, seperated by spaces or commas" )
    ]
    def _init( self, *args, **kwargs ):
        self.pdb_id_list = map( 
            lambda x: x[0:4], 
            re.split("[\s,]+", self.pdb_id.strip() ) 
        )
        self.pdb_file_list = map(
            lambda x: self.outpath( "%s.pdb" % x ),
            self.pdb_id_list
        )    
        self.output_files += self.pdb_file_list
    def func( self ):
        for pdb_id, pdb_file in zip( self.pdb_id_list, self.pdb_file_list ):
            pdb_download( pdb_id, pdb_file )




def pdb_assembly( pdb_id ):
    try:
        url = PDBE_ASSEMBLY.format( pdb_id=pdb_id )
        return urllib2.urlopen( url ).read()
    except urllib2.HTTPError:
        raise Exception( "Error downloading pdb assembly '%s'" % pdb_id )

class PdbAssembly( PyTool ):
    args = [
        _( "pdb_id", type="text", 
            help="pdb id" )
    ]
    out = [
        _( "assembly_file", file="assembly.pdb" )
    ]
    def func( self ):
        with open( self.assembly_file, "w" ) as fp:
            fp.write( pdb_assembly( self.pdb_id ) )




def pdb_split( pdb_file, output_dir, backbone_only=False, 
               resno_ignore=False, max_models=False, zfill=False ):
    """ author: Johanna Tiemann
        author: Alexander Rose
        This function puts pdb-models into their own file.
    """ 
    backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
    bb_tag = "_bb" if backbone_only else ""
    model_no = 0
    if max_models and not zfill:
        zfill = len( str(max_models) )
    with open( pdb_file, "r" ) as fp:
        for line in fp:
            if line[0:5]=='MODEL':
                model_no += 1
                if max_models and model_no>max_models:
                    break
                file_name = "%s%s.pdb" % ( str(model_no).zfill( zfill ), bb_tag )
                file_path = os.path.join( output_dir, file_name )
                with open( file_path, 'w') as fp_out:
                    for line in fp:
                        if line[0:4]!='ATOM':
                            if line[0:6]=='ENDMDL':
                                fp_out.write( 'END' )
                                break
                            continue
                        if backbone_only and line[12:16] not in backbone:
                            continue
                        if resno_ignore:
                            if try_int( line[22:26] ) in resno_ignore:
                                continue
                        fp_out.write( line )
                        
class PdbSplit( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "backbone_only", type="checkbox", default=False ),
        _( "max_models", type="slider", range=[0, 100], default=0 ),
        _( "resno_ignore", type="text", default="" ),
        _( "zfill", type="slider", range=[0, 8], default=0 )
    ]
    def _init( self, *args, **kwargs ):
        if self.resno_ignore:
            if isinstance( self.resno_ignore, basestring ):
                self.resno_ignore = map( int, self.resno_ignore.split(",") )
    def func( self ):
        pdb_split( 
            self.pdb_file, self.output_dir, 
            backbone_only=self.backbone_only,
            max_models=self.max_models,
            resno_ignore=self.resno_ignore,
            zfill=self.zfill
        )
class LoopDelete (PyTool):
    args = [
    _( "pdb_file", type="file", ext="pdb" ),
    _( "chain", type="string"),
    _( "res1", type="int"),
    _( "res2", type="int"),
    ]
    out = [
    _( "noloop_pdb_file", file="noloop.pdb" )
    ]
    def func( self ):
        npdb = numpdb.NumPdb( self.pdb_file, {
            "phi_psi": False,
            "sstruc": False,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False,
            "info": False
        })
        test =npdb.sele(chain=self.chain,resno=[ self.res1, self.res2], invert=True)
        
        npdb.write( self.noloop_pdb_file,sele=test )
        
class PdbEdit( PyTool ):
    """ Edits a pdb file. Manipulations are done in this order:
        center, shift, box.
    """
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
        _( "center", type="checkbox", default=False ),
        _( "shift", type="float", nargs=3, default=None,
            metavar=("X", "Y", "Z") ),
        _( "box", type="float", nargs=6, default=None,
            help="Two sets of 3d coordinates; first one corner, "
                "then the extent at which the opposing corner is. "
                "Outputs only the atoms within the box.",
            metavar=("X", "Y", "Z", "EX", "EY", "EZ") ),
    ]
    out = [
        _( "edited_pdb_file", file="edited.pdb" )
    ]
    def func( self ):
        npdb = numpdb.NumPdb( self.pdb_file, {
            "phi_psi": False,
            "sstruc": False,
            "backbone_only": False,
            "protein_only": False,
            "detect_incomplete": False,
            "configuration": False,
            "info": False
        })
        sele = None

        if self.center:
            npdb['xyz'] -= npdb['xyz'].mean( axis=0 )

        if self.shift:
            shift = np.array( self.shift[0:3] )
            npdb['xyz'] += shift

        if self.box:
            corner1 = np.array( self.box[0:3] )
            corner2 = corner1 + np.array( self.box[3:6] )
            sele = (
                ( npdb['xyz']>corner1 ) & 
                ( npdb['xyz']<corner2 )
            ).all( axis=1 )

        npdb.copy( sele=sele ).write( self.edited_pdb_file )


class PdbSuperpose( PyTool, ProviMixin ):
    args = [
        _( "pdb_file1", type="file", ext="pdb" ),
        _( "pdb_file2", type="file", ext="pdb" ),
        _( "sele1", type="sele" ),
        _( "sele2", type="sele" ),
        _( "subset|ss", type="text", default="CA" ),
        _( "rmsd_cutoff|co", type="float", default=1.0 )
    ]
    out = [
        _( "superposed_file", file="superposed.pdb" )
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "superpose.provi"
    def func( self ):
        npdb1 = numpdb.NumPdb( self.pdb_file1 )
        npdb2 = numpdb.NumPdb( self.pdb_file2 )
        numpdb.superpose( 
            npdb1, npdb2, self.sele1, self.sele2, 
            subset=self.subset, inplace=True,
            rmsd_cutoff=self.rmsd_cutoff, max_cycles=100
        )
        npdb1.write( self.superposed_file )
    def _post_exec( self ):
        self._make_provi_file(
            pdb_file1=self.relpath( self.pdb_file1 ),
            pdb_file2=self.relpath( self.pdb_file2 ),
            superposed_file=self.relpath( self.superposed_file ),
        )



class PdbInfo( PyTool ):
    args = [
        _( "pdb_input", type="str" )
    ]
    out = [
        _( "info_file", file="pdb_info.json" )
    ]
    def _init( self, *args, **kwargs ):
        if len(self.pdb_input)==4:
            self.pdb_download = PdbDownload(
                self.pdb_input,
                **copy_dict( kwargs, run=False,
                    output_dir=self.subdir("download") )
            )
            self.pdb_file = self.pdb_download.pdb_file_list[0]
        else:
            self.pdb_download = None
            self.pdb_file = self.abspath( self.pdb_input )
    def _pre_exec( self ):
        if self.pdb_download:
            self.pdb_download()
    def func( self ):
        self.info = numpdb.NumPdb( self.pdb_file, features={
            "phi_psi": False, 
            "info": True,
            "backbone_only": True
        })._info
        with open( self.info_file, "w" ) as fp:
            json.dump( self.info, fp, indent=4 )
    def get_info( self ):
        if not hasattr( self, "info" ):
            with open( self.info_file, "r" ) as fp:
                self.info = json.load( fp )
        return self.info


def numpdb_test( pdb_file ):
    with Timer("read/parse pdb plain"):
        numpdb.NumPdb( pdb_file, features={
            "phi_psi": False, 
            "sstruc": False, 
            "backbone_only": True,
            "detect_incomplete": False
        })
    with Timer("read/parse pdb incomplete"):
        numpdb.NumPdb( pdb_file, features={"phi_psi": False, "sstruc": False, "backbone_only": False, "detect_incomplete": False} )
    with Timer("read/parse pdb phi/psi"):
        numpdb.NumPdb( pdb_file, features={"sstruc": False} )
    with Timer("read/parse pdb"):
        npdb = numpdb.NumPdb( pdb_file )
    # with Timer("resno iter2"):
    #     for numa_list in npdb.iter_resno2( 6, chain="A", resno=[1,22] ):
    #         print [ numa["resno"][0] for numa in numa_list ]
    # with Timer("dist"):
    #     print npdb.dist( {"chain":"A"}, {"chain":"B"} )
    with Timer("access phi/psi"):
        print np.nansum( npdb['phi'] ), np.nansum( npdb['psi'] )
    with Timer("access phi/psi, altloc"):
        print np.nansum( npdb.get('phi', altloc=[" ", "A"] ) )
    with Timer("sstruc iter"):
        for numa in npdb.iter_sstruc():
            pass
    with Timer("sequence"):
        print npdb.sequence(chain="A", resno=[1,100])
    
    npdb = numpdb.NumPdb( pdb_file, features={
        "phi_psi": False, 
        "sstruc": False, 
        "backbone_only": False,
        "detect_incomplete": False
    })
    with Timer("plain resno iter"):
        for numa in npdb._iter_resno():
            pass
    with Timer("write pdb"):
        npdb.write( "test.pdb" )

    npdb = numpdb.NumPdb( pdb_file, features={
        "phi_psi": False, 
        "sstruc": False, 
        "backbone_only": False,
        "detect_incomplete": True
    })
    with Timer("resno iter new"):
        for numa in npdb.iter_resno():
            pass


class NumpdbTest( PyTool ):
    args = [
        _( "pdb_file", type="file", ext="pdb" )
    ]
    def func( self ):
        numpdb_test( self.pdb_file )




def rcsb_search( data ):
    req = urllib2.Request( PDB_SEARCH_URL, data=data )
    result = urllib2.urlopen(req).read()
    return result.split() if result else []

def rna_list( max_res=3.5 ):
    """ author: Johanna Tiemann, Alexander Rose
        query rcsb for pdb files containing rna
    """
    date = '1955-01-01'
    today = str( datetime.datetime.now().strftime("%Y-%m-%d") )
    pdb_list = []
    searchstr_res = [
        'Resolution',
        'refine.ls_d_res_high.min>0.0</refine.ls_d_res_high.min',
        'refine.ls_d_res_high.max>'+str(max_res)+'</refine.ls_d_res_high.max'
    ]
    searchstr_nmr = [
        'ExpType',
        ('mvStructure.expMethod.value>'
            'SOLUTION NMR'
        '</mvStructure.expMethod.value'),
        ('mvStructure.hasExperimentalData.value>'
            'Y'
        '</mvStructure.hasExperimentalData.value')
    ]
    rna_query_temp = string.Template(
        open( os.path.join( TMPL_DIR, "rna_query.tpl" ) ).read()
    )
    for q in [ searchstr_res, searchstr_nmr ]:
        query_text = rna_query_temp.substitute(
            ins1=q[0], ins2=q[1], ins3=q[2],
            date=date, today=today
        )
        pdb_list += rcsb_search( query_text )
    list_record = ListRecord(
        "rna_query", PDB_SEARCH_URL,
        None, None, pdb_list
    )
    return list_record


#############ROTAMERE-START#############


# the rotamere lib needs the 'remaining_atoms' data
ROTAMERE_LIB_PATH = os.path.join(
    BASEKIT_DIR, "basekit","data", "bio", "bbind02.May.lib.json"
)
with open( ROTAMERE_LIB_PATH, "r" ) as fp:
    ROTAMERE_LIB = json.load( fp )



def get_rotno( resname ):
    # get number of available rotameres for that resname
    rota=ROTAMERE_LIB.get('dihedral_angles')
    rnmbr=len(rota.get(resname))
    return rnmbr

def get_rotamere( resname, no ):
    # return dihedral_angle, dihedral_atoms, remaining_atoms
    rota=ROTAMERE_LIB.get('dihedral_angles')
    rotanrdi=rota.get(resname)[no]
    dihi=ROTAMERE_LIB.get('dihedral_atoms')
    dihiat=dihi.get(resname)
    remat=ROTAMERE_LIB.get('remaining_atoms')
    remaining_atoms=remat.get(resname)
    return rotanrdi, dihiat, remaining_atoms
def rmake_rotamere(   npdb, sele, no ):
    dihedral_angle, dihedral_atoms, remaining_atoms =  get_rotamere( sele["resname"], no )
    #print sele
    am1=npdb.sele(resname=sele["resname"],resno=sele["resno"],chain=sele["chain"])
    print "sle",  am1
    am=npdb.copy(sele=am1)
    print am['xyz']
    for chi_index in range( 0,len( dihedral_angle )-1 ):
        test=am.sele(atomname=huhu)#(atomname=dihedral_atoms[chi_index])
        print test
        hurz=am.copy(sele=test)
        print dihedral_atoms[chi_index]
        print hurz['xyz']
def make_rotamere( npdb, sele, no ):
    print  sele, no
    dihedral_angle, dihedral_atoms, remaining_atoms =  get_rotamere( sele["resname"], no )
    # geht ueber alle chi-angle (im Beispiel CYS nur einmal, da nur ein chi-angle)
    coords = npdb.get( 'xyz', **sele )
    print 'alle coordinaten:', coords
    for chi_index in range( 0,len( dihedral_angle )-1 ):#len( dihedral_angle )-1 ):
        print 'chi angle_NR:',chi_index
        
        
        atom_pos = npdb.get( 'atomname', **sele )
        atom_no = npdb.get( 'atomno', **sele )
        print 'alle coordinaten:', coords
        print 'alle coordinaten:', atom_pos
        print 'alle coordinaten:', atom_no
        
        # TODO: ptr = points to rotate --> oder alle rein?
        #       waere dann coords
        
        # TODO: wtr = what should be rotated: defined by 0,1...
        #       1 defines the coords, who are holded
        
        
        #coords of remaining atoms:
        coords_remat=np.empty( [len( remaining_atoms[ chi_index ] ), 3] )
        coords_diheat=np.empty( [4, 3] )
        i=0;j=0;
        for index, coord_atom in enumerate( atom_pos ):
            for diheat in dihedral_atoms[ chi_index ]:
                if diheat.split() == coord_atom.split():
                    coords_diheat[ i ] = coords[ index ]
                    i+=1
                    break;
            for remat in remaining_atoms[ chi_index ]:
                if remat == coord_atom.split()[0]:
                    coords_remat[j]=coords[ index ]
                    j+=1
                    break;
            
        #print 'Koordinaten der remaining atoms:', coords_remat
        #print 'Koordinaten der dihedral atoms:', coords_diheat
        
        
        curr_dihedral=numpdb.dihedral(coords_diheat[0],coords_diheat[1],coords_diheat[2],coords_diheat[3])
        #print dihedral(13.112 13.937 )
        av_angle = dihedral_angle[ chi_index+1 ]
        #print'#####', dihedral([15.036,  11.747,  18.715],[13.564,  11.573,  18.836],[12.933,  12.737,  19.580],[13.206,  12.785,  21.061])
        
        
        
        #rotation auf 180 gesetz um es besser visualisieren zu koennen
        
        rotation =av_angle-curr_dihedral
        if rotation <0:
            rotation =360+rotation
        print 'Current dihedral angle', curr_dihedral
        print 'rotation angles', rotation
        print 'Dihedral angle of ROTAMERE_LIB:', av_angle
        
        newpoints = np.empty( [len( remaining_atoms[ chi_index ] ), 3] )
        for index, remat_co in enumerate(coords_remat):
            oldpoint=np.array(remat_co)
            wtr=coords_diheat[2]-coords_diheat[1]
            shift=coords_diheat[2]*-1
            oldpoint=oldpoint+shift
            rotmat = rmatrixu( wtr, np.deg2rad( rotation ) )#3.14159
            newpoint = ( np.dot( rotmat, oldpoint.T ) ).T
            newpoint=newpoint-shift
            newpoints[index] = newpoint
        newpdb=np.copy(coords)
        j=0
        for coord_atom in atom_pos:
            for index, remat in enumerate(remaining_atoms[ chi_index ]):
                if remat == coord_atom.split()[0]:
                    newpdb[j]=newpoints[index]
            j+=1

        #print 'newpdb ', newpdb
        coords = np.copy(newpdb)
        print 'alle coordinaten2:', coords
    #update npdb
    all_coords = npdb['xyz']
    for index, at_num in enumerate(atom_no):
        all_coords[at_num-1] = newpdb[index]
    npdb['xyz'] = all_coords
    
    
    print 'alle coordinaten:', npdb.get( 'atomname', **sele )
    #print 'updated coords', npdb['xyz']
    
    
    #
    
    
    
    
    return npdb#rotamer_pdb
    #
        
        #
        #kjljl=jlkj
        #for index, remat in enumerate(remaining_atoms[ chi_index ]):
        #    if remat == atoms['atomname'].split()[0]:
        #        newpdb[atind]['xyz']=newpoints[index]
        #
        #
        #newpdb=np.copy(coords)
        #atom_pos = npdb.get( 'atomname', **sele )
        #for atind, atoms in enumerate(newpdb):
        #    for index, remat in enumerate(remaining_atoms[ chi_index ]):
        #        if remat == atoms['atomname'].split()[0]:
        #            print newpdb['xyz']
        #            newpdb[atind]['xyz']=newpoints[index]
        #        
        
        #newpdb=np.copy(coords)
        #j=0
        #for coord_atom in atom_pos:
        #    for index, remat in enumerate(remaining_atoms[ chi_index ]):
        #        if remat == coord_atom.split()[0]:
        #            newpdb[j]=newpoints[index]
        #    j+=1
        #npdb = np.copy(newpdb)
        #print 'newpdb ', newpdb
        #for elem in npdb:
        ##    print elem
        #
        #coords['xyz']= coords
        #
        #for elem in npdb:
        #    print elem
        #
        #
        #
        
        #update coords
        #coords=newpoints
       # npdb.copy( sele=sele ).write( self.edited_pdb_file )
    
    pass
    #makes the rotation and gives the new coords back
    #rotmat = rmatrixu( wtr, np.deg2rad( rotangl ) )
    #newpoints = ( np.dot( rotmat, ptr.T ) ).T
    
    #TODO: put newpoints into npdb
    

    
    
    # use rmatrixu from utils.math, see test/test_utils_math.py
    # for an example how to use the rotation matrix
    # - for all chis (mind the order)
    #   - calc current dihedral
    #   - calc how much the dihedral needs to be rotated
    #   - rotate...
    # return npdb with changed rotamere
    #return npdbr
# put the "make all rotameres" function in pdb.py

def make_all_rotameres(pdbfile, chain1, resno1,zfill):
    print pdbfile, chain1, resno1
    rotamere_dict={}
    npdb = numpdb.NumPdb( pdbfile)
    sele={"resno": resno1, "chain": chain1}
    resname1=npdb.get('resname', **sele)[0]
    print 'resname',resname1
    sele={ "resno": resno1, "chain": chain1, "resname": resname1 }
    print sele
    print npdb.get( 'xyz', **sele )
    no = get_rotno( sele["resname"] )
    test =npdb.sele(chain=chain1,resno=resno1)
    for i in range(0, no):
        rotamere = make_rotamere( npdb, sele, i )
        rotamere.copy(sele=test).write("%s_%s_%s.%s" % (resname1,resno1, str(i+1).zfill(zfill) ,'pdb'))
        #npdb.copy(rotamer)
        #rotamere_dict[ (pdbfile, chain1, resno1, i) ] = rotamere
        print rotamere
        #npdb.write(test.pdb)
    return rotamere_dict, test
    
class MakeAllRotameres( PyTool ):
    args = [
        _( "pdb_input", type="file", ext="pdb" ),
        _( "chain|ch", type="text" ),
        _( "resno|r", type="int" ),
        _( "zfill", type="slider", range=[0, 8], default=0 )
    ]
    out = [
        _( "rotamere_file", file="{pdb_input.stem}.pdb" )
    ]
    def func( self, *args, **kwargs ):
        print self.chain, self.resno
        make_all_rotameres( self.pdb_input, self.chain, self.resno ,self.zfill)
    #def _post_exec( self ):
    #    self.rotamere_record
    #    for index, elem in enumerate(self.rotamere_record):
    #        self.rotamere_record[elem].write(self.rotamere_file+'.'+self.chain+'.'+str(self.resno)+'.'+str(index)+'.pdb', sele=self.test)



#############ROTAMERE-END#############

class RnaList( PyTool ):
    args = [
        _( "compare_list|cp", type="file", ext="json", default='' ),
        _( "max_resolution|mr", type="float", default=3.5 )
    ]
    out = [
        _( "current_list", file="current.json" ),
        _( "new_list", file="new.json" ),
        _( "old_list", file="old.json" )
    ]
    def func( self ):
        self.rna_record = rna_list( self.max_resolution )
    def _post_exec( self ):
        ListIO( self.current_list ).write( self.rna_record )
        if self.compare_list:
            new_record, old_record = list_compare(
                self.rna_record, ListIO( self.compare_list ).read()
            )
            ListIO( self.new_list ).write( new_record )
            ListIO( self.old_list ).write( old_record )




class ListCompare( PyTool ):
    # TODO make it work for an arbitrary number of lists
    #   by doing all pairwise comparisions
    args = [
        _( "list1", type="file", ext="json" ),
        _( "list2", type="file", ext="json" ),
        _( "list_name|ln", type="text", default=None )
    ]
    out = [
        _( "list1_only", file="list1_only.json" ),
        _( "list2_only", file="list2_only.json" )
    ]
    def func( self ):
        list1_only, list2_only = list_compare(
            ListIO( self.list1 ).read(),
            ListIO( self.list2 ).read(),
            name=self.list_name
        )
        ListIO( self.list1_only ).write( list1_only )
        ListIO( self.list2_only ).write( list2_only )


class ListJoin( PyTool ):
    args = [
        _( "list", type="file", ext="json", nargs="+" ),
        _( "list_name|ln", type="text", default=None )
    ]
    out = [
        _( "joined_list", file="joined_list.json" )
    ]
    def func( self ):
        joined_list = list_join(
            *map( lambda x: ListIO( x ).read(), self.list ),
            name=self.list_name
        )
        ListIO( self.joined_list ).write( joined_list )




