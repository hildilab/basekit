from __future__ import with_statement
from __future__ import division



import os

from utils import copy_dict
from utils.tool import CmdTool
from utils.numpdb import NumPdb, numsele

import provi_prep as provi



DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]

try:
    LINKIT_DIR = os.environ["LINKIT_DIR"]
except:
    print "LINKIT_DIR environment variable not defined"
    LINKIT_DIR = ""
LINKIT_CMD = os.path.join( LINKIT_DIR, "Link_It_dos2.exe" )




class Linker( CmdTool ):
    args = [
        { "name": "pdb_file", "type": "file", "ext": "pdb" },
        { "name": "res1", "type": "text" },
        { "name": "res2", "type": "text" },
        { "name": "seq", "type": "text" }
    ]
    def _init( self, pdb_file, res1, res2, seq, **kwargs ):
        self.pdb_file = os.path.abspath( pdb_file )
        stem = os.path.splitext( os.path.split( self.pdb_file )[-1] )[0]
        self.bin_file = os.path.join( self.output_dir, "%s_linker.bin" % stem )
        self.txt_file = os.path.join( self.output_dir, "%s_linker.txt" % stem )
        self.kos_file = os.path.join( self.output_dir, "%s_kos.txt" % stem )
        self.res1 = res1
        self.res2 = res2
        self.seq = seq
        self.cmd = [ "wine", LINKIT_CMD, self.kos_file, self.bin_file, "t" ]
        self.output_files = [ 
            self.bin_file, self.txt_file
        ]
    def _make_kos( self ):
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
                print sele, atomname, coords
                fp.write( "%s\n" % "\n".join(map( str, coords[0] ) ) )
            fp.write( "%s\n" % self.seq )
            for sele in ( sele1, sele2 ):
                fp.write( "%s\t%s\n" % ( sele.get("chain", " "), sele["resno"] ) )
    def _pre_exec( self ):
        self._make_kos()



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