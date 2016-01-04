from __future__ import with_statement
from __future__ import division


import os
import re
import json
import urllib2
import xml.etree.ElementTree

from poster.encode import multipart_encode
from poster.streaminghttp import register_openers

import numpy as np

from utils import memoize_m, try_float
from utils.math import norm
from utils.tool import _, _dir_init, PyTool, ProviMixin
from utils.listing import ListRecord, ListIO
import utils.numpdb as numpdb
from utils.numpdb import NumPdb, numsele

DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "tmdet" )
PDBTM_LOCAL_PATH = os.environ.get("PDBTM_LOCAL_PATH", "")
TMDET_URL = "http://tmdet.enzim.hu/"
PDBTM_ALPHAHELICAL_URL = "http://pdbtm_data.enzim.hu/pdbtmalpha"
PDBTM_ALPHAHELICAL_PATH = os.environ.get("PDBTM_ALPHAHELICAL_PATH", "")


def pdbtm_download():
    try:
        url = PDBTM_ALPHAHELICAL_URL
        return urllib2.urlopen( url ).read()
    except urllib2.HTTPError:
        raise Exception("MPstruc url error")


def pdbtm( pdbid, db_path=PDBTM_LOCAL_PATH ):
    """queries a local PDBTM data base given a four character pdb id"""
    with open( db_path ) as db:
        # <pdbtm xmlns="http://pdbtm.enzim.hu" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://pdbtm.enzim.hu/pdbtm.xsd pdbtm.xsd" ID="1ap9" TMP="yes">
        start_re = re.compile( '<pdbtm .* ID="%s".*>' % pdbid )
        inside_entry = False
        result = ''
        for line in db:
            if not inside_entry:
                if start_re.search(line):
                    if line.find('TMP="yes"') >= 0:
                        inside_entry = True
                        result += line
                    else:
                        # found an entry which is not an membrane protein
                        return False
                        break
            else:
                result += line
                if line.find('</pdbtm>') >= 0:
                    inside_entry = False
                    break;
        if result:
            return '<?xml version="1.0" encoding="iso-8859-1"?>\n' + result
        else:
            raise Exception('No results in tmdet db found.')




def tmdet( pdb_file ):
    """queries the TMDET webservice"""

    # Register the streaming http handlers with urllib2
    register_openers()

    # headers contains the necessary Content-Type and Content-Length
    # datagen is a generator object that yields the encoded parameters
    datagen, headers = multipart_encode({
        "pdbfile": open( pdb_file ),
        "go": "cgi"
    })

    # Create the Request object
    request = urllib2.Request(TMDET_URL + "index.php", datagen, headers)
    response = urllib2.urlopen(request).read()

    no_tmp = re.search('is not a transmembrane protein', response)
    if no_tmp:
        return False

    # HREF="tmp/phpd0dDXC.xml"
    m = re.search('HREF="(tmp/php.*\.xml)"', response)

    if m:
        responseXml = urllib2.urlopen( TMDET_URL + m.group(1) )
        return responseXml.read()
    else:
        raise Exception('No results in tmdet webservice response found.')


class PdbtmDownload( PyTool ):
    """A tool to download the Pdbtm database"""
    out = [
        _( "pdbtm_xml", file="pdbtm.xml" )
    ]
    def func( self ):
        with open( self.pdbtm_xml, "w" ) as fp:
            fp.write( pdbtm_download() )



class TmdetMixin( object ):
    def _post_exec( self ):
        return
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



class Tmdet( TmdetMixin, PyTool, ProviMixin ):
    """A tool to access the TMDET server"""
    args = [
        _( "pdb_file", type="file" ),
    ]
    out = [
        _( "tmdet_file", file="{pdb_file.stem}.xml" ),
        #_( "mplane_file", file="{pdb_file.stem}.mplane" ),
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "tmdet.provi"
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        tmdet_xml = tmdet( self.pdb_file )
        if not tmdet_xml:
            tmdet_xml=" "
        with open( self.tmdet_file, "w" ) as fp:
            fp.write( tmdet_xml )
        with open(self.tmdet_file, 'r') as fp:
            text = fp.read()
            if text != " ":
                O = []
                N = []
                s = text
                s = s.replace("\n", " ")
                result = re.findall('<CHAIN(.*?)TYPE="alpha"', s)
                if result:
                    chain = re.findall('CHAINID="(.*?)"', result[0])
                    regiontext = re.findall('TYPE="alpha"(.*?)</CHAIN>',s)
                    regions = re.findall('<REGION (.*?)type="H"/>', regiontext[0])
                    found = []
                    for reg in regions:
                        regio = reg.split("REGION")[-1]
                        pdb_beg = re.findall('pdb_beg="(.*?)"', regio)
                        pdb_end = re.findall('pdb_end="(.*?)"', regio)
                        found.append([ int(pdb_beg[0])-1, int(pdb_end[0])+1 ])
                    for index, elem in enumerate(found):
                        N.append(elem[index%2])
                        O.append(elem[(index+1)%2])
                    npdb = NumPdb( self.pdb_file, features={
                        "phi_psi": False, "sstruc": False, "backbone_only": True
                    })
                    sele = {}
                    
                    
                    with open(self.output_dir + self.pdb_file.split('.')[0].split('/')[-1]+'_tmdet.ply', 'w') as fp2:
                        fp2.write(  "ply\n"+\
                                    "format ascii 1.0\n"+\
                                    "element vertex 8\n"+\
                                    "property float x\n"+\
                                    "property float y\n"+\
                                    "property float z\n"+\
                                    "property uchar red\n"+\
                                    "property uchar green\n"+\
                                    "property uchar blue\n"+\
                                    "element face 4\n"+\
                                    "property list uchar int vertex_indices\n"+\
                                    "property uchar red\n"+\
                                    "property uchar green\n"+\
                                    "property uchar blue\n"+\
                                    "end_header\n")
                        oldx = []
                        oldy = []
                        oldz = []
                        for val in [N, O]:
                            alln = []
                            if val == N:
                                sele["atomname"] = "C"
                            elif val == O:
                                sele["atomname"] = "N"
                            else:
                                print 'error'
                            for elem in val:
                                sele["resno"] = elem
                                coords = npdb.get( 'xyz', **sele )
                                alln.append(coords[0].tolist())
                            x = [j[0] for j in alln]
                            y = [j[1] for j in alln]
                            z = [j[2] for j in alln]
                            if val == N:
                                oldx = x
                                oldy = y
                                oldz = z
                            if val == O:
                                #x = oldx
                                y = oldy
                                z = oldz
                            print np.std(x)
                            fp2.write( str(np.mean(x))+" "+str( (np.mean(y)+np.std(y))*1.0)+" "+str( (np.mean(z)+np.std(z))*1.0)+" 255\t0\t0\n")
                            fp2.write( str(np.mean(x))+" "+str( (np.mean(y)-np.std(y))*1.0)+" "+str( (np.mean(z)+np.std(z))*1.0)+" 255\t0\t0\n")
                            fp2.write( str(np.mean(x))+" "+str( (np.mean(y)+np.std(y))*1.0)+" "+str( (np.mean(z)-np.std(z))*1.0)+" 255\t0\t0\n")
                            fp2.write( str(np.mean(x))+" "+str( (np.mean(y)-np.std(y))*1.0)+" "+str( (np.mean(z)-np.std(z))*1.0)+" 255\t0\t0\n")
                        fp2.write("3 0 1 2 255\t0\t0\n")
                        # fp2.write("3 0 1 3 255\t0\t0\n")
                        # fp2.write("3 0 2 3 255\t0\t0\n")
                        fp2.write("3 1 2 3 255\t0\t0\n")
                        fp2.write("3 4 5 6 255\t0\t0\n")
                        # fp2.write("3 4 5 7 255\t0\t0\n")
                        # fp2.write("3 4 6 7 255\t0\t0\n")
                        fp2.write("3 5 6 7 255\t0\t0\n")





def pdbtm_tree( xml_file=None ):
    if not xml_file:
        xml_file = PDBTM_ALPHAHELICAL_PATH
    try:
        tree = xml.etree.ElementTree.parse( xml_file )
    except IOError:
        raise Exception( "PDBTM xml file not found: '%s'" % xml_file )
    return tree.getroot()



class PdbtmDb( object ):
    def __init__( self, xml_file=None ):
        self.tree = pdbtm_tree( xml_file=xml_file )
    def find( self, pdb_id ):
        pdb_id = pdb_id.lower()
        for protein in self.tree:
            if protein.attrib["ID"].lower()==pdb_id:
                return protein
    def list( self ):
        pdbid_list = []
        for protein in self.tree:
            pdbid_list.append( protein.attrib["ID"] )
        return pdbid_list


class Pdbtm( TmdetMixin, PyTool, ProviMixin ):
    pass


def pdbtm_list( xml_file=None ):
    pm = PdbtmDb( xml_file=xml_file )
    list_record = ListRecord(
        "pdbtm", PDBTM_ALPHAHELICAL_PATH,
        None, None, pm.list()
    )
    return list_record


class PdbtmList( PyTool ):
    args = [

    ]
    out = [
        _( "list_file", file="pdbtm_list.json" )
    ]
    def func( self ):
        list_record = pdbtm_list()
        ListIO( self.list_file ).write( list_record )



def pdbtm_info(  ):
    pass

class PdbtmInfo( PyTool ):
    """A tool to get infos from the PDBTM database"""
    args = [
        _( "pdb_id", type="str" ),
    ]
    out = [
        _( "info_file", file="pdbtm_info_{pdb_id}.json" )
    ]
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        self.info = pdbtm_info( self.pdb_id )
        with open( self.info_file, "w" ) as fp:
            json.dump( self.info, fp, indent=4 )
    def get_info( self ):
        if not hasattr( self, "info" ):
            with open( self.info_file, "r" ) as fp:
                self.info = json.read( fp )
        return self.info

