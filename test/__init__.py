
import os


DIR = os.path.split( os.path.abspath(__file__) )[0]
PARENT_DIR = os.path.split( DIR )[0]
DATA_DIR = os.path.join( PARENT_DIR, "data", "test" )

TMP_DIR = os.path.join( DIR, "tmp" )



# cd ./test/
# python -m unittest discover


