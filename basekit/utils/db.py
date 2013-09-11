import os
import sqlite3



def get_pdb_subdirs( directory ):
    subdirs = []
    for subdir in os.listdir( directory ):
        subdir_path = os.path.join( directory, subdir )
        if os.path.isdir( subdir_path ):
            subdirs.append( subdir_path )
    subdirs.sort()
    return subdirs


def get_pdb_subdir_files( subdir, pattern=None ):
    l = []
    for fname in os.listdir( subdir ):
        if not pattern or fname.endswith( pattern ):
            fpath = os.path.join( subdir, fname )
            l.append( fpath )
    return l


def get_pdb_files( directory, pattern=None ):
    subdirs = get_pdb_subdirs( directory )
    pdb_files = []
    for subdir in subdirs:
        pdb_files += get_pdb_subdir_files( subdir, pattern=pattern )
    return pdb_files




def table_existance( c, name ):
    print name, len(name)
    c.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name=?", [name])
    return c.fetchone()[0]


def create_table( db_path, schema, name, data, overwrite=False ):
    with sqlite3.connect( db_path, isolation_level="EXCLUSIVE" ) as conn:
        c = conn.cursor()

        table_exist = table_existance( c, name )

        if table_exist:
            if overwrite:
                c.execute("DROP TABLE %s" % name)
            else:
                conn.commit()
                print "table '%s' exists and no overwrite cmd given, exiting" % name
                return

        # create table
        c.execute( schema )

        # populate table
        val = ",".join( ['?']*len( data[0] ) )
        c.executemany('INSERT INTO %s VALUES (%s)' % (name, val), data)

        conn.commit()


