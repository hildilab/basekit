from __future__ import division

import itertools

import numpy as np
np.seterr( all="raise" )


# https://github.com/pycogent/pycogent/blob/master/cogent/struct/dihedral.py



def vec_norm( v ):
    if v.shape == (3,):
        return v/vec_mag(v)
    else:
        mag = vec_mag(v)
        v2 = np.copy(v)
        v2[:,0] /= mag
        v2[:,1] /= mag
        v2[:,2] /= mag
        return v2


def vec_mag( v ):
    if v.shape == (3,):
        return np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    else:
        return np.sqrt( v[:,0]**2 + v[:,1]**2 + v[:,2]**2 )


def vec_dot( v1, v2 ):
    if v1.shape == (3,) and v2.shape == (3,):
        return np.dot( v1, v2 )
    else:
        return np.sum( v1*v2, axis=1 )


def vec_angle( v1, v2 ):
    if v1.shape == (3,) and v2.shape == (3,):
        dot = np.dot(v1, v2)
    elif v1.shape == (3,) and v2.shape != (3,):
        dot = np.dot(v1, v2.T)
    elif v1.shape != (3,) and v2.shape == (3,):
        dot = np.dot(v1, v2)
    else:
        n = max(len(v1), len(v2))
        dot = np.array([ np.dot(v1[i], v2[i]) for i in range(n) ])
    ang = np.arccos( dot / (vec_mag(v1)*vec_mag(v2)) )
    ang = ang*180 / np.pi
    return ang


def vec_project2plane( v, pn ):
    c1 = np.cross(v, pn)
    c2 = np.cross(c1, pn)
    return vec_norm(c2)


def vec_dihedral( v1, v2, v3, v4 ):
    v12 = v2-v1
    v23 = v3-v2
    v34 = v4-v3

    n1 = vec_norm( np.cross( v12, v23 ) )
    n2 = vec_norm( np.cross( v23, v34 ) )

    torsion = vec_angle( n1, n2 )
    if vec_dot( n1, v34 ) < 0:
        torsion *= -1
    return torsion



def norm( v ):
    return v/np.sqrt( np.dot( v, v ) )

def mag( v ):
    return np.sqrt( np.dot( v, v ) )
    
def angle( v1, v2 ):
    dot = np.dot( v1, v2 )
    mag = np.sqrt( np.dot( v1, v1 ) * np.dot( v2, v2 ) )
    ang = np.arccos( dot / mag )
    return ang*180 / np.pi

def dihedral( v1, v2, v3, v4 ):
    v12 = v2-v1
    v23 = v3-v2
    v34 = v4-v3

    c1 = np.cross( v12, v23 )
    c2 = np.cross( v23, v34 )

    torsion = angle( c1, c2 )

    if np.dot( c1, v34 ) < 0:
        torsion *= -1
    return torsion




def lsq(y):
    # y = mx + c
    x = np.arange( len(y) )
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq( A, y )[0]
    return [ m*x[0]+c, m*x[-1]+c ]

def axis( coords ):
    return np.array([ lsq(coords[...,i]) for i in range(3) ]).T

def rmsd( coords1, coords2 ):
    return np.sqrt( 
        np.sum( np.power( coords1-coords2, 2 ) ) / coords1.shape[0]
    )







class Superposition( object ):
    def __init__( self, src_coords, dst_coords ):
        self.src_coords = src_coords
        self.dst_coords = dst_coords
        self.n = src_coords.shape[0]
        self.rmsd = 0.0
        self.quaternion = None
        self.rotmat = None
        self._superpose()
    def _superpose( self ):
        self.src_center = self.src_coords.mean( axis=0 )
        self.dst_center = self.dst_coords.mean( axis=0 )
        v0 = ( self.src_coords.copy() - self.src_center ).T
        v1 = ( self.dst_coords.copy() - self.dst_center ).T

        # SVD of covar matrix
        u, s, vh = np.linalg.svd( np.dot( v1, v0.T ) )
        # rotation matrix from SVD orthonormal bases
        R = np.dot( u, vh )
        if np.linalg.det( R ) < 0.0:
            # R not a right handed system
            R -= np.outer( u[:, 2], vh[2, :]*2.0 )
            s[-1] *= -1.0
        # homogeneous transformation matrix
        M = np.identity(4)
        M[:3, :3] = R        

        # translation
        M[:3, 3] = self.dst_center
        T = np.identity(4)
        T[:3, 3] = -self.src_center
        M = np.dot(M, T)

        # rotation matrix
        self.rotmat = M[0:3,0:3]

        # rmsd
        E0 = np.sum(np.sum( v0*v0 )) + np.sum(np.sum( v1*v1 ))
        msd = ( E0 - 2.0*sum(s) ) / v0.shape[1]
        self.rmsd = np.sqrt( max([ msd, 0.0 ]) )
    def transform( self, coords, inplace=False ):
        if not inplace:
            position = coords.copy()
        else:
            position = coords
        position -= self.src_center
        position = ( np.dot( self.rotmat, position.T ) ).T
        position += self.dst_center
        return position





