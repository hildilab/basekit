import os
import sys
import unittest

DIR = os.path.split( os.path.abspath( __file__ ) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )
TMP_DIR = os.path.join( DIR, "tmp" )

sys.path.append( os.path.join( PARENT_DIR, "data", "ball" ) )

# Examples: http://ball-trac.bioinf.uni-sb.de/wiki/CodeLibrary
from BALL import *



def data( file_name ):
    return os.path.join( DATA_DIR, file_name )

def tmp( *dir_name ):
    return os.path.join( TMP_DIR, "BALL", *dir_name )


class BALLTestCase( unittest.TestCase ):
    def setUp( self ):
        pass
    def test_check( self ):
        # read the PDB-file into a BALL::System
        f = PDBFile( data( "1CRN.pdb" ) )
        S = System()
        f.read(S) 

        # now we open a fragment database
        fdb = FragmentDB( "" )

        # and normalize the atom names, i.e. we convert different
        # naming standards to the PDB naming scheme - just in case!
        S.apply( fdb.normalize_names )

        # now we add any missing hydrogens to the residues
        # the data on the hydrogen positions stems from the
        # fragment database. However the hydrogen positions
        # created in this way are only good estimates
        S.apply( fdb.add_hydrogens )

        # now we create the bonds between the atoms (PDB files hardly
        # ever contain a complete set of CONECT records)
        S.apply( fdb.build_bonds )

        # check the first protein
        protein = S.getProtein( 0 )

        if( protein != None ):                      
            self.assertEqual(
                Peptides.GetSequence( protein ),
                "TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN"
            )
            self.assertEqual( protein.countChains(), 1 )
            self.assertEqual( chains( protein )[0].countResidues(), 46 )

