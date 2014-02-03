import os

import moderna as modeRNA
import numpy as np
from utils.tool import _, _dir_init, PyTool
import os



DIR, PARENT_DIR, TMPL_DIR = _dir_init( __file__, "moderna_tools" )



def load_templates( pdbfile, chain='A' ):
    """
    Loads a template structure from a PDB file. Produces a Template object that can be saved in a variable.
    Each template in ModeRNA has only one chain. By default, the chain with id 'A' is loaded. Another chain id can be specified optionally
    Arguments:
        path+filename of a PDB structure file
        chain id (by default 'A')
        """
    
    return modeRNA.load_template( pdbfile, chain )

def load_models( pdbfile, chain='A', data_type='file', template=None, alignment=None ):
    """
    Loads a structure model that has been built previously, or any PDB struture that is then to be modified. Produces a RnaModel object that can be saved in a variable. Each model in ModeRNA contains only one chain. Multi-chain structures can be modeled by using more than one Template/Alignment/RNAModel at a time.
    By default, RNAModels are created by reading files. They can also be created from BioPython PDB objects (precisely Bio.PDB.Structure.Structure and Bio.PDB.Chain.Chain objects)
    Arguments:	
        path+filename of a PDB structure file; or a Structure or Chain object from BioPython.
        chain id (by default 'A')
        data type ('file' by default; if set to 'structure' or 'chain', Bio.PDB objects are read)
        Template object to be used for this model (optional)
        Alignment object to be used for this model (optional)
        """
    return modeRNA.load_model( pdbfile, chain, template, alignment )

def finds_modifications( structure ):
    """
    Reports all modified nucleotides that occur in a structure loaded or created with ModeRNA (models, templates).
    This command returns all modified residues as a Python dictionary.
        {'1':<res>, '6': <res>, '6a': <res>, ...}
    Arguments:
        RnaModel or Template object
        """
    return modeRNA.find_modifications( structure )

def examine_structures( structure, logfile=False ):
    """
    Checks whether the given structure has any features that may cause any problems during the modeling process. The user needs to give a structure object, and optionally a name of the file the report is written to.
    Arguments:
        Stucture object
        name of logfile (optional)
        """
    return modeRNA.examine_structure( structure, logfile )

def clean_structures( structure, write_structures=True ):
    """
    Eliminates features that may cause problems during modeling from a template or model structure:
        water molecules, ions, amino acids, and unidentified residues are deleted.
        old atom names are replaced (C1* becomes C1', O1P becomes OP1).
        missing phosphate groups are added (also adds the OP3 atom for these)
        reports whether the chain is continuous.
    In case some feature cannot be fixed (e.g. chain discontinuity) this is written to the logfile. It is recommended to include such features in the alignment ('.' characters for strange residues, and '_' for backbone breaks).
    Arguments:
        Stucture object (RnaModel or Template)
        True/False - whether structure should be written to a PDB file (optional)
        """
    return modeRNA.clean_structure( structure, write_structure=write_structures )

def removes_all_modifications( model ):
    """
    Removes all base modifications from a given model. The nucleotides are transformed into standard bases from which the modifications originated.
    Arguments:
        RnaModel object
        """
    return modeRNA.remove_all_modifications( model )

def add_pdbs( original_file, pdb_file, moderna_file, chain ):
    def write_to_file(read_file, write_file, chain):
        with open(write_file, "a") as fp:
                with open(read_file, "r") as fp2:
                    for line in fp2:
                        if line.startswith("ATOM") and line[21]==chain:
                            fp.write(line)
    try:
        if os.stat(pdb_file).st_size > 0:
            write_to_file(pdb_file, moderna_file, chain)
            with open(pdb_file, "w"):
                pass
        else:
            write_to_file(original_file, moderna_file, chain)
            print chain
    except OSError:
        write_to_file(original_file, moderna_file, chain)
    


def finish_pdbs( moderna_file ):
    with open(moderna_file, "a") as fp:
        fp.write("TER")

class Moderna( PyTool ):
    """
    MODERNA is a selection of tools to manipulate the pdb.
        """
    args = [
        _( "pdb_input", type="file" ),
        _( "logfile|l", type="bool", default=False,
            help="Default: False" ),
        _( "write_structure|ws", type="bool", default=True,
            help="Writes after cleaning the pdb-structure out. Default: False" ),
        
        _( "tool", type="str",
            options=[   "finds_modifications", "examine_structures",
                        "clean_structures", "removes_all_modifications"
                    ], default="finds_modifications",
            help="What should be done. options: finds_modifications,"+
                    "examine_structures, clean_structures, "+
                    "removes_all_modifications. Default: finds_modifications" )
    ]
    out = [
        _( "moderna_file_temp", file="fixed_structure.pdb" ),
        _( "moderna_file", file="moderna.pdb" ),
        _( "log_file", file="moderna.log" )
    ]
    tmpl_dir = TMPL_DIR

    def func( self ):
        sele = []
        with open(self.pdb_input, "r") as fp:
            for line in fp:
                if line.startswith("ATOM"):
                    sele.append(line[21])
        sele = np.array(sele)
        chains = np.unique(sele)
        if self.tool=='finds_modifications':
            for chain in chains:
                pdbfile = load_templates( self.pdb_input, chain )
                finds_modifications( pdbfile )
                add_pdbs( self.pdb_input, self.moderna_file_temp, self.moderna_file, chain )
            finish_pdbs( self.moderna_file )
        if self.tool=='examine_structures':
            for chain in chains:
                pdbfile = load_templates( self.pdb_input, chain )
                examine_structures( pdbfile, self.logfile )
                add_pdbs( self.pdb_input, self.moderna_file_temp, self.moderna_file, chain )
            finish_pdbs( self.moderna_file )
        if self.tool=='clean_structures':
            for chain in chains:
                pdbfile = load_templates( self.pdb_input, chain )
                clean_structures( pdbfile, self.write_structure )
                add_pdbs( self.pdb_input, self.moderna_file_temp, self.moderna_file, chain )
            finish_pdbs( self.moderna_file )
        if self.tool=='removes_all_modifications':
            for chain in chains:
                pdbfile = load_models( self.pdb_input, chain )
                removes_all_modifications( pdbfile )
                add_pdbs( self.pdb_input, self.moderna_file_temp, self.moderna_file, chain )
            finish_pdbs( self.moderna_file )
        
        
        