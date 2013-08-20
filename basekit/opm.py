from __future__ import with_statement
from __future__ import division


import os
import re
import json
import urllib2

import numpy as np


from utils import memoize_m, try_float
from utils.math import norm
from utils.tool import _, _dir_init, PyTool, ProviMixin



DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "opm" )
OPM_URL = "http://opm.phar.umich.edu/"
OPM_PDB_URL = OPM_URL + "pdb/{pdb_id:}.pdb"
OPM_LIST_URL = OPM_URL + "superfamilies_dl.php?class={class_id:}"
OPM_SEARCH_URL = OPM_URL + "protein.php?search={pdb_id:}"
OPM_LOCAL_PATH = os.environ.get("OPM_LOCAL_PATH", "")


def opm_list( class_id ):
    try:
        url = OPM_LIST_URL.format( class_id=class_id )
        page = urllib2.urlopen( url ).read()
    except urllib2.HTTPError:
        raise Exception("Opm_list url error")
    pdb_ids = re.findall( r"([0-9a-zA-Z]{4})<br />", page )
    pdb_ids = [ x.upper() for x in pdb_ids ]
    pdb_ids = list( set( pdb_ids ) ) # make unique
    pdb_ids.sort()
    return pdb_ids



def opm_info( pdb_id ):
    try:
        url = OPM_SEARCH_URL.format( pdb_id=pdb_id )
        page = urllib2.urlopen( url ).read()
    except urllib2.HTTPError:
        raise Exception("Opm_info url error")

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
        "OPMType": opm_type[0].split(" ", 1)[1],
        "Class": opm_class[0].split(" ", 1)[1],
        "OPMSuperfamily": opm_superfamily[0].split(" ", 1)[1],
        "OPMFamily": opm_family[0].split(" ", 1)[1],
        "OPMSpecies": opm_species[0].split(" ", 1)[1],
        "OPMLocalization": opm_localization[0],
        "OPMRelated": related_ids, 
        "OPMDeltaG": try_float( delta_g[0] )
    }




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


def ppm( pdb_file ):
    raise NotImplementedError


class OpmMixin( object ):
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
        # TODO works only for two planes
        with open( self.opm_file ) as fp:
            coords={'O':[],'N':[]}
            for line in fp:
                if line[0:6]=="HETATM" and line[17:20]=="DUM":
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
                raise Exception( "could not find plane coordinates" )
        return np.array([ coords["N"], coords["O"] ])
    def make_mplane_file( self ):
        mp = self.get_planes().tolist()
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
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "opm.provi"
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        with open( self.opm_file, "w" ) as fp:
            fp.write( opm( self.pdb_id ) )



class Ppm( OpmMixin, PyTool, ProviMixin ):
    """A tool to query the PPM webservice"""
    args = [
        _( "pdb_file", type="file", ext="pdb" ),
    ]
    out = [
        _( "opm_file", file="{pdb_file.stem}_opm.pdb" ),
        _( "mplane_file", file="{pdb_file.stem}.mplane" ),
        _( "processed_file", file="{pdb_file.stem}_proc.pdb" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "opm.provi"
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        ppm( self.pdb_file )



class OpmList( PyTool ):
    """A tool to access the OPM database"""
    args = [
        _( "class_id", type="str" ),
    ]
    out = [
        _( "list_file", file="opm_list_class_{class_id}.txt" )
    ]
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        pdb_list = opm_list( self.class_id )
        with open( self.list_file, "w" ) as fp:
            fp.write( "\n".join( pdb_list ) )


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
        info = opm_info( self.pdb_id )
        with open( self.info_file, "w" ) as fp:
            json.dump( info, fp, indent=4 )

