from __future__ import with_statement
from __future__ import division


import os
import json
import urllib2
import xml.etree.ElementTree

from utils.tool import _, _dir_init, PyTool
from utils.list import ListRecord, ListIO



DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "mpstruc" )
MPSTRUC_URL = "http://blanco.biomol.uci.edu/mpstruc/"
MPSTRUC_ALPHAHELICAL_URL = MPSTRUC_URL + "listAll/mpstrucAlphaHlxTblXml"
MPSTRUC_ALPHAHELICAL_PATH = os.environ.get("MPSTRUC_ALPHAHELICAL_PATH", "")



def mpstruc_download():
    try:
        url = MPSTRUC_ALPHAHELICAL_URL
        return urllib2.urlopen( url ).read()
    except urllib2.HTTPError:
        raise Exception("MPstruc url error")


def mpstruc_tree( xml_file=None ):
    if not xml_file:
        xml_file = MPSTRUC_ALPHAHELICAL_PATH
    try:
        tree = xml.etree.ElementTree.parse( xml_file )
    except IOError:
        raise Exception( "MPstruc xml file not found: '%s'" % xml_file )
    return tree.getroot()



class MpstrucDb( object ):
    def __init__( self, xml_file=None ):
        self.tree = mpstruc_tree( xml_file=xml_file )
    def find( self, pdb_id ):
        pdb_id = pdb_id.lower()
        for protein in self.tree.findall(".//protein"):
            if protein.find("pdbCode").text.lower()==pdb_id:
                return protein
            rel = protein.find("relatedPdbEntries")
            if rel!=None:
                for x in rel.findall("pdbCode"):
                    if x.text.lower()==pdb_id:
                        return protein
        for protein in self.tree.findall(".//memberProtein"):
            if protein.find("pdbCode").text.lower()==pdb_id:
                return protein
            rel = protein.find("relatedPdbEntries")
            if rel!=None:
                for x in rel.findall("pdbCode"):
                    if x.text.lower()==pdb_id:
                        return protein
    def _group( self, pdb_id, sub=False ):
        key = ".//subgroup" if sub else ".//group"
        pdb_id = pdb_id.lower()
        for group in self.tree.findall( key ):
            for protein in group.findall(".//protein"):
                if protein.find("pdbCode").text.lower()==pdb_id:
                    return group.find("name").text
                rel = protein.find("relatedPdbEntries")
                if rel!=None:
                    for x in rel.findall("pdbCode"):
                        if x.text.lower()==pdb_id:
                            return group.find("name").text
            for protein in group.findall(".//memberProtein"):
                if protein.find("pdbCode").text.lower()==pdb_id:
                    return group.find("name").text
    def group( self, pdb_id ):
        return self._group( pdb_id )
    def subgroup( self, pdb_id ):
        return self._group( pdb_id, sub=True )
    def master( self, pdb_id ):
        protein = self.find( pdb_id )
        master_tag = protein.find("masterProteinPdbCode")
        if master_tag!=None:
            protein = self.find( master_tag.text )
        return protein
    def info( self, pdb_id ):
        protein = self.find( pdb_id )
        if protein!=None:
            return {
                "name": protein.find("name").text,
                "species": protein.find("species").text,
                "subgroup": self.subgroup( pdb_id ),
                "group": self.group( pdb_id )
            }
        else:
            return None
    def str( self, protein ):
        for x in protein:
            print "%s: %s" % ( x.tag, x.text )
    def list( self ):
        pdbid_list = []
        for protein in self.tree.findall(".//protein"):
            pdbid_list.append( protein.find("pdbCode").text )
            rel = protein.find("relatedPdbEntries")
            if rel!=None:
                for x in rel.findall("pdbCode"):
                    pdbid_list.append( x.text )
        for protein in self.tree.findall(".//memberProtein"):
            pdbid_list.append( protein.find("pdbCode").text )
            rel = protein.find("relatedPdbEntries")
            if rel!=None:
                for x in rel.findall("pdbCode"):
                    pdbid_list.append( x.text )
        return [ x.upper() for x in pdbid_list ]


def mpstruc_info( pdb_id, xml_file=None ):
    pdb_id = pdb_id.lower()
    mp = MpstrucDb( xml_file=xml_file )
    return mp.info( pdb_id )



class MpstrucDownload( PyTool ):
    """A tool to download the MPstruc database"""
    out = [
        _( "mpstruc_xml", file="mpstruc.xml" )
    ]
    def func( self ):
        with open( self.mpstruc_xml, "w" ) as fp:
            fp.write( mpstruc_download() )


class MpstrucInfo( PyTool ):
    """A tool to get infos from the MPstruc database"""
    args = [
        _( "pdb_id", type="str" ),
        _( "mpstruc_xml|mx", type="file", ext="xml", default="" ),
    ]
    out = [
        _( "info_file", file="mpstruc_info_{pdb_id}.json" )
    ]
    def _init( self, *args, **kwargs ):
        pass
    def func( self ):
        self.info = mpstruc_info( self.pdb_id, xml_file=self.mpstruc_xml )
        with open( self.info_file, "w" ) as fp:
            json.dump( self.info, fp, indent=4 )
    def get_info( self ):
        if not hasattr( self, "info" ):
            with open( self.info_file, "r" ) as fp:
                self.info = json.load( fp )
        return self.info



def mpstruc_list( xml_file=None ):
    mp = MpstrucDb( xml_file=xml_file )
    list_record = ListRecord(
        "mpstruc", mp.tree.attrib["url"],
        None, mp.tree.attrib["lastDatabaseEditDate"], mp.list()
    )
    return list_record


class MpstrucList( PyTool ):
    args = [
        
    ]
    out = [
        _( "list_file", file="mpstruc_list.json" )
    ]
    def func( self ):
        list_record = mpstruc_list()
        ListIO( self.list_file ).write( list_record )



