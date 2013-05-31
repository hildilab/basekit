from __future__ import with_statement
from __future__ import division



import os
import itertools
import operator
import json

import utils.path
from utils import copy_dict, iter_stride
from utils.tool import CmdTool, PyTool, ProviMixin
from utils.numpdb import NumPdb, numsele

import provi_prep as provi
from spider import LoopCrosscorrel



DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
TMPL_DIR = os.path.join( PARENT_DIR, "data", "linker" )

try:
    LINKIT_DIR = os.environ["LINKIT_DIR"]
except:
    print "LINKIT_DIR environment variable not defined"
    LINKIT_DIR = ""
LINKIT_CMD = os.path.join( LINKIT_DIR, "Link_It_dos2.exe" )



class LinkerTest( PyTool ):
    args = [
        { "name": "linker_txt", "type": "file", "ext": "txt" }
    ]
    def _init( self, linker_txt, **kwargs ):
        self.linker_txt = self.abspath( linker_txt )
        self.linker_json = self.outpath( "%s.json" % utils.path.stem( linker_txt, "json" ) )
    def func( self ):
        self._make_linker_json( compact=True )
    def _make_linker_json( compact=False ):
        linker_dict = {}
        with open( self.linker_txt, "r" ) as fp:
            x = fp.next()
            n = fp.next()
            for i, d in enumerate( iter_stride( fp, 3 ), start=1 ):
                linker_dict[ i ] = [ float(d[0]), float(d[1]), d[2].strip() ]
        with open( self.linker_json, "w" ) as fp:
            if compact:
                json.dump( linker_dict, fp, separators=(',',':') )
            else:
                json.dump( linker_dict, fp, indent=4 )


class LinkIt( CmdTool, ProviMixin ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "res1", "type": "text" },
        { "name": "res2", "type": "text" },
        { "name": "seq", "type": "text" }
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it.provi"
    def _init( self, pdb_file, res1, res2, seq, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        stem = utils.path.stem( self.pdb_file )
        self.bin_file = os.path.join( self.output_dir, "%s_linker.bin" % stem )
        self.txt_file = os.path.join( self.output_dir, "%s_linker.txt" % stem )
        self.pdb_linker_file = os.path.join( self.output_dir, "%s_linker.pdb" % stem )
        self.pdb_linker_file2 = os.path.join( self.output_dir, "%s_linker2.pdb" % stem )
        self.pdb_linker_file3 = os.path.join( self.output_dir, "%s_linker3.pdb" % stem )
        self.kos_file = os.path.join( self.output_dir, "%s_kos.txt" % stem )
        self.res1 = res1
        self.res2 = res2
        self.seq = seq
        self.json_file = self.outpath( "%s.json" % utils.path.stem( self.txt_file ) )
        self.cmd = [ "wine", LINKIT_CMD, self.kos_file, self.bin_file, "t" ]
        self.output_files = [ 
            self.bin_file, self.txt_file, self.pdb_linker_file, self.pdb_linker_file2
        ]
    def _pre_exec( self ):
        self._make_kos_file()
    def _post_exec( self ):
        self._fix_linker_pdb( self.pdb_linker_file2 )
        self._fix_linker_pdb( self.pdb_linker_file3, atoms_only=True )
        self._make_linker_json( compact=True )
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            pdb_linker_file3=self.relpath( self.pdb_linker_file3 ),
            json_file=self.relpath( self.json_file )
        )
    def _make_kos_file( self ):
        npdb = NumPdb( self.pdb_file, features={ 
            "phi_psi": False, "sstruc": False, "backbone_only": True
        })
        sele1 = numsele( self.res1 )
        sele2 = numsele( self.res2 )
        with open( self.kos_file, "w" ) as fp:
            d = [ (sele1, " CA "), (sele1, " C  "), (sele2, " N  "), (sele2, " CA ") ]
            for sele, atomname in d:
                sele["atomname"] = atomname
                coords = npdb.get( 'xyz', **sele )
                fp.write( "%s\n" % "\n".join(map( str, coords[0] ) ) )
            fp.write( "%s\n" % self.seq )
            for sele in ( sele1, sele2 ):
                fp.write( "%s\t%s\n" % ( sele.get("chain", " "), sele["resno"] ) )
    def _fix_linker_pdb( self, output_file, atoms_only=False, stems=True ):
        backbone = ( ' N  ',' C  ', ' CA ',' O  ' )
        with open( self.pdb_linker_file, "r" ) as fp:
            with open( output_file, "w" ) as fp_out:
                for i, line in enumerate( fp ):
                    if line.startswith("MODEL"):
                        atom_i = 1
                    if line.startswith("ATOM"):
                        line = line = line[0:6] + ( "% 5i" % atom_i ) + line[11:]
                        if line[22]=="X":
                            if not stems:
                                continue
                            if line[12:16] not in backbone:
                                continue
                            tag = "1000 " if line[24]==" " else "2000 "
                            line = line[0:17] + "GLY" + line[20:22] + tag + line[27:]
                        atom_i += 1
                        fp_out.write( line )
                        continue
                    if not atoms_only:
                        fp_out.write( line )
    def _make_linker_json( self, compact=False ):
        linker_dict = {}
        with open( self.txt_file, "r" ) as fp:
            x = fp.next()
            n = fp.next()
            for i, d in enumerate( iter_stride( fp, 3 ), start=1 ):
                linker_dict[ i ] = [ float(d[0]), float(d[1]), d[2].strip() ]
        with open( self.json_file, "w" ) as fp:
            if compact:
                json.dump( linker_dict, fp, separators=(',',':') )
            else:
                json.dump( linker_dict, fp, indent=4 )

    

class LinkItDensity( PyTool, ProviMixin ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "mrc_file", "type": "file", "ext": "mrc" },
        { "name": "res1", "type": "text" },
        { "name": "res2", "type": "text" },
        { "name": "seq", "type": "text" },
        { "name": "pixelsize", "type": "slider", "range": [1, 10], "fixed": True },
        { "name": "resolution", "type": "slider", "range": [1, 10], "fixed": True },
        { "name": "max_loops", "type": "slider", "range": [0, 200], "default_value": 100 }
    ]
    tmpl_dir = TMPL_DIR
    provi_tmpl = "link_it_density.provi"
    def _init( self, pdb_file, mrc_file, res1, res2, seq, 
               pixelsize, resolution, max_loops=100, **kwargs ):
        self.pdb_file = self.abspath( pdb_file )
        self.mrc_file = self.abspath( mrc_file )
        self.link_it = LinkIt( 
            self.pdb_file, res1, res2, seq,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("link_it") )
        )
        self.loop_correl = LoopCrosscorrel(
            self.mrc_file, self.pdb_file, self.link_it.pdb_linker_file2, 
            res1, res2, len(seq), pixelsize, resolution, max_loops=max_loops,
            **copy_dict( kwargs, run=False, output_dir=self.subdir("loop_correl") )
        )
        self.output_files = list( itertools.chain(
            self.link_it.output_files, 
            self.loop_correl.output_files,
        ))
    def func( self ):
        self.link_it()
        self.loop_correl()
    def _post_exec( self ):
        self._make_provi_file(
            pdb_file=self.relpath( self.pdb_file ),
            mrc_file=self.relpath( self.mrc_file ),
            box_mrc_file=self.relpath( self.loop_correl.spider_reconvert.mrc_file ),
            pdb_linker_file3=self.relpath( self.link_it.pdb_linker_file3 ),
            link_it_json_file=self.relpath( self.link_it.json_file )
        )
    def _make_correl_json( self, compact=False ):
        linker_dict = {}
        with open( self.txt_file, "r" ) as fp:
            x = fp.next()
            n = fp.next()
            for i, d in enumerate( iter_stride( fp, 3 ), start=1 ):
                linker_dict[ i ] = [ float(d[0]), float(d[1]), d[2].strip() ]
        with open( self.json_file, "w" ) as fp:
            if compact:
                json.dump( linker_dict, fp, separators=(',',':') )
            else:
                json.dump( linker_dict, fp, indent=4 )



# 1. File analog dem angehaengten schreiben (KOS.txt):
# 12 Zeilen Koordinaten, gewuenschte Sequenz, gewuenschte
# Aminosaeure-Numerierung (erste/letzte)
# 2. dos.exe aufrufen: 'Link_it_dos.exe KOS.txt MyLinker.bin t'
# (t steht fuer txt-Variante, also nicht binaer)
# 3. 'Link_it_dos.exe rechnet, schreibt das Ausgabefile, und hinterher ein
# File 'Ready.bin', deshalb:
# Warten, bis ready.bin aufgetaucht ist, danach linkerfile lesen

# KO_TxtFile := ExtractFilePath(ParamStr(0))+'KOS.txt';
# Linker_Antwort := ExtractFilePath(ParamStr(0))+'MyLinker.bin';
# Ready_File := ExtractFilePath(ParamStr(0))+'Ready.bin';
# Link_It_Dos_Path := ExtractFilePath(ParamStr(0));
# Link_it_Web_Ini_File:=ExtractFilePath(ParamStr(0))+'link_it_web.ini';
# Linker_Antwort := ThePath+'MyLinker.bin';


# Aufruf:='Link_it_dos.exe '+KO_TxtFile+' '+Linker_Antwort+' t';
# GotLinker:=False;
# WinExecAndWait_32(Aufruf,Link_It_Dos_Path, SW_SHOW, True);

# while not(FileExists(Ready_File)) do Application.HandleMessage;
# while not(FileExists(Linker_Antwort)) do Application.ProcessMessages;




# Link_it.ini: Pfade zu db setzen



# kos.txt:
# 12 Zeilen Koordinaten
# 1-3 CA 1.AS
# 4-6 C   1.AS
# 7-9 N    2.AS
# 10-12 CA 2.AS
# Sequenz des Linkers
# Stem1  Stem2 (mit Kette)


# wine Link_It_dos2.exe KOS.txt mylinker.bin t


# mylinker.bin => output
# t => output additionally as text