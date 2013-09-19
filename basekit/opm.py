from __future__ import with_statement
from __future__ import division


import os
import re
import json
import urllib2
import collections

from poster.encode import multipart_encode
from poster.streaminghttp import register_openers

import numpy as np

from utils import memoize_m, try_float, copy_dict, get_index
from utils.math import norm
from utils.tool import _, _dir_init, PyTool, ProviMixin
from utils.list import ListRecord, ListIO
from utils.numpdb import NumPdb

from pdb import PdbAssembly


DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "opm" )
OPM_URL = "http://opm.phar.umich.edu/"
OPM_PDB_URL = OPM_URL + "pdb/{pdb_id:}.pdb"
OPM_LIST_URL = OPM_URL + "superfamilies_dl.php?class={class_id:}"
OPM_SEARCH_URL = OPM_URL + "protein.php?search={pdb_id:}"
OPM_LOCAL_PATH = os.environ.get("OPM_LOCAL_PATH", "")
PPM_URL = "http://sunshine.phar.umich.edu/"


def opm_list( class_id ):
    try:
        url = OPM_LIST_URL.format( class_id=class_id )
        page = urllib2.urlopen( url ).read()
    except urllib2.HTTPError:
        raise Exception("Opm_list url error")
    pdb_ids = re.findall( r"([0-9a-zA-Z]{4})<br />", page )
    list_record = ListRecord(
        "opm superfamily %s" % class_id, url,
        None, None, pdb_ids
    )
    return list_record


def _parse_opm_info( page ):

    # check if there were no matches
    no_matches = re.findall(
        r'<h2>Search Results for ".*"</h2>No matches', page
    )
    if no_matches:
        return None

    # check if this page only points to a representative structure
    rep = re.findall(
        r'Representative structure\(s\) of this protein: <br /> '
        r'<a href="protein\.php\?pdbid=([0-9a-zA-Z]{4})">', page
    )
    if rep:
        return { "representative": rep[0].upper() }

    opm_type = re.findall( 
        r'<li><i>Type:</i> <a.*>(.*)</a>', page )
    opm_class = re.findall( 
        r'<li><i>Class:</i> <a.*>(.*)</a>', page )
    opm_superfamily = re.findall( 
        r'<li><i>Superfamily:</i> <a[^<]*>([^<]*)</a>', page )
    opm_family = re.findall( 
        r'<li><i>Family:</i> <a[^<]*>([^<]*)</a>', page )
    opm_species = re.findall( 
        r'<li><i>Species:</i> <i><a.*>(.*)</a></i>', page )
    opm_localization = re.findall( 
        r'<li><i>Localization:</i> <a.*>(.*)</a>', page )

    related_ids = re.findall( r'"\?extrapdb=([0-9a-zA-Z]{4})"', page )
    related_ids = [ x.upper() for x in related_ids ]
    related_ids.sort()

    delta_g = re.findall( r'([-+]?[0-9]*\.?[0-9]+) kcal/mol', page )

    return {
        "type": opm_type[0].split(" ", 1)[1],
        "class": opm_class[0].split(" ", 1)[1],
        "superfamily": opm_superfamily[0].split(" ", 1)[1],
        "family": opm_family[0].split(" ", 1)[1],
        "species": opm_species[0].strip(),
        "localization": opm_localization[0],
        "related_ids": related_ids, 
        "delta_g": try_float( get_index( delta_g, 0 ) )
    }


def opm_info( pdb_id ):
    try:
        url = OPM_SEARCH_URL.format( pdb_id=pdb_id )
        page = urllib2.urlopen( url ).read()
    except urllib2.HTTPError:
        raise Exception("Opm_info url error")

    with open( "test_opm_info.html", "w" ) as fp:
        fp.write( page )

    info = _parse_opm_info( page )
    if info:
        rep = info.get( "representative", None )
        # get info from representative entry
        if rep:
            info_rep = opm_info( rep )
            info = info_rep
            if pdb_id in info["related_ids"]:
                info["related_ids"].remove( pdb_id )
            info["related_ids"].append( rep )
            info["representative"] = rep

    return info




def opm( pdb_id ):
    pdb_id = pdb_id.lower()
    try:
        path = os.path.join( OPM_LOCAL_PATH, "%s.pdb" % pdb_id )
        with open( path, "r" ) as fp:
            return fp.read()
    except IOError:
        try:
            url = OPM_PDB_URL.format( pdb_id=pdb_id )
            return urllib2.urlopen( url ).read()
        except urllib2.HTTPError:
            raise Exception("Opm url error")


def _parse_ppm( page ):
    msg_list = [
        r'returned an error: (.*)$',
        r'(Too many residues)'
    ]
    for msg in msg_list:
        error = re.findall( msg, page )
        if error:
            with open( "ppm_error.txt", "w" ) as fp:
                fp.write( error[0] )
            raise Exception( error[0] )
    pdb_url = PPM_URL + re.findall( 
        r'href="\./(pdb_upload/.*out\.pdb)"', page
    )[0]
    delta_g = re.findall( r'([-+]?[0-9]*\.?[0-9]+) kcal/mol', page )[0]
    info_dict = {
        "delta_g": try_float( delta_g )
    }
    return pdb_url, info_dict


def ppm( pdb_file, topology="in", hetero=True ):
    """queries the PPM webservice"""
    # Register the streaming http handlers with urllib2
    register_openers()
    datagen, headers = multipart_encode({
        "inout": "in" if topology else "out",
        "yesno": "yes" if hetero else "no",
        "userfile": open( pdb_file ),
        "submit": "Submit"
    })
    request = urllib2.Request( PPM_URL + "upload_file.php", datagen, headers )
    page = urllib2.urlopen( request ).read()
    
    with open( "test_ppm.html", "w" ) as fp:
        fp.write( page )
    # with open( "test_ppm.html", "r" ) as fp:
    #     page = fp.read()
    pdb_url, info_dict = _parse_ppm( page )
    pdb_file = urllib2.urlopen( pdb_url ).read()

    info_dict.update({
        "hetero": hetero,
        "topology": topology    
    })
    return pdb_file, info_dict


def parse_planes( opm_pdb_file ):
    with open( opm_pdb_file ) as fp:
        coords = collections.OrderedDict([ ('N', []), ('O',[]) ])
        for line in fp:
            if line[0:6] in ( "ATOM  ", "HETATM" ) and line[17:20]=="DUM":
                atm = line[13:14]
                c = np.array( map( float, [ 
                    line[30:38], line[38:46], line[46:54] 
                ]))
                if len(coords[ atm ])==2:
                    vn1 = norm( coords[ atm ][1] - coords[ atm ][0] )
                    vn2 = norm( c - coords[ atm ][0] )
                    if not np.allclose( vn1, vn2 ):
                        coords[ atm ].append(c)
                elif len(coords[ atm ])<2:
                    coords[ atm ].append(c)
            if sum( map( len, coords.values() ) )==6:
                break
        else:
            if len( coords["N"] )==3:
                del coords["O"]
            elif len( coords["O"] )==3:
                del coords["N"]
            else:
                raise Exception( "could not find plane coordinates" )
    return np.array( coords.values() )


class OpmMixin( object ):
    tmpl_dir = TMPL_DIR
    provi_tmpl = "opm.provi"
    def _post_exec( self ):
        with open( self.processed_file, "w" ) as fp:
            with open( self.opm_file, "r" ) as fp_opm:
                for line in fp_opm:
                    if line[17:20]!="DUM" and line[0:6]!="REMARK":
                        fp.write( line )
        self.make_mplane_file()
        self._make_provi_file(
            pdb_file=self.relpath( self.processed_file ),
            mplane_file=self.relpath( self.mplane_file )
        )
    @memoize_m
    def get_planes( self ):
        with open( self.mplane_file, "r" ) as fp:
            return json.load( fp )
    def make_mplane_file( self ):
        mp = parse_planes( self.opm_file ).tolist()
        if mp:
            with open( self.mplane_file, "w" ) as fp:
                json.dump( mp, fp )



class Opm( OpmMixin, PyTool, ProviMixin ):
    """A tool to access the OPM database pdb files"""
    args = [
        _( "pdb_id", type="str" ),
    ]
    out = [
        _( "opm_file", file="{pdb_id}_opm.pdb" ),
        _( "mplane_file", file="{pdb_id}.mplane" ),
        _( "processed_file", file="{pdb_id}_proc.pdb" ),
        _( "info_file", file="{pdb_id}_info.json", optional=True ),
    ]
    def func( self ):
        with open( self.opm_file, "w" ) as fp:
            fp.write( opm( self.pdb_id ) )



class PpmMixin( OpmMixin ):
    def ppm( self, pdb_input_file ):
        pdb_file, info_dict = ppm( pdb_input_file )
        with open( self.opm_file, "w" ) as fp:
            fp.write( pdb_file )
        with open( self.info_file, "w" ) as fp:
            json.dump( info_dict, fp, indent=4 )


class Ppm( PpmMixin, PyTool, ProviMixin ):
    """A tool to query the PPM webservice with a pdb file"""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
    ]
    out = [
        _( "opm_file", file="{pdb_file.stem}_ppm.pdb" ),
        _( "mplane_file", file="{pdb_file.stem}.mplane" ),
        _( "processed_file", file="{pdb_file.stem}_proc.pdb" ),
        _( "info_file", file="{pdb_file.stem}_info.json" ),
    ]
    def func( self ):
        self.ppm( self.pdb_file )
    


class Ppm2( PpmMixin, PyTool, ProviMixin ):
    """A tool to query the PPM webservice with a pdb id"""
    args = [
        _( "pdb_id", type="str" ),
    ]
    out = [
        _( "opm_file", file="{pdb_id}_ppm.pdb" ),
        _( "mplane_file", file="{pdb_id}.mplane" ),
        _( "assembly_file", file="{pdb_id}_asm.pdb" ),
        _( "processed_file", file="{pdb_id}_proc.pdb" ),
        _( "info_file", file="{pdb_id}_info.json" ),
    ]
    def _init( self, *args, **kwargs ):
        self.pdb_assembly = PdbAssembly(
            self.pdb_id,
            **copy_dict( kwargs, run=False, 
                output_dir=self.subdir("assembly") )
        )
    def _pre_exec( self ):
        self.pdb_assembly()
        # extract first model
        npdb = NumPdb( 
            self.pdb_assembly.assembly_file, 
            features={
                "phi_psi": False, 
                "info": False,
            }
        )
        npdb.write2( self.assembly_file )
    def func( self ):
        self.ppm( self.assembly_file )


class OpmList( PyTool ):
    args = [
        _( "class_id", type="str" ),
    ]
    out = [
        _( "list_file", file="opm_list_class_{class_id}.json" )
    ]
    def func( self ):
        list_record = opm_list( self.class_id )
        ListIO( self.list_file ).write( list_record )


class OpmInfo( PyTool ):
    """A tool to get infos from the OPM database"""
    args = [
        _( "pdb_id", type="str" ),
    ]
    out = [
        _( "info_file", file="opm_info_{pdb_id}.json" )
    ]
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        self.info = opm_info( self.pdb_id )
        with open( self.info_file, "w" ) as fp:
            json.dump( self.info, fp, indent=4 )
    def get_info( self ):
        if not hasattr( self, "info" ):
            with open( self.info_file, "r" ) as fp:
                self.info = json.load( fp )
        return self.info

