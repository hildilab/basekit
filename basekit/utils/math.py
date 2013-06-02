from __future__ import division

import numpy as np


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





