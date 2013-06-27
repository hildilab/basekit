#include "Python.h"
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <math.h>

// crossProduct3, dotProduct3, dihedral mostly from 
// https://github.com/SimTk/msmbuilder/blob/master/src/ext/dihedral/dihedral_wrap.c

inline void crossProduct3(double a[], const double b[], const double c[]) {
    //Calculate the cross product between length-three vectors b and c, storing
    //the result in a
    (a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2];
    (a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0];
    (a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];
}

inline double dotProduct3(const double b[], const double c[]) {
    //Calculate the dot product between length-three vectors b and c
    return b[0] * c[0] + b[1] * c[1] + b[2] * c[2];
}

double dihedral(const double *x0, const double *x1, const double *x2,
                const double *x3) {
    //Calculate the signed dihedral angle between four points.
    //Result in radians
    //x0, x1, x2, x3 should be length three arrays
    int i;
    double b1[3], b2[3], b3[3], c1[3], c2[3];
    double arg1, arg2, b2_norm;
    
    for (i = 0; i < 3; i++) {
        b1[i] = x1[i] - x0[i];
        b2[i] = x2[i] - x1[i];
        b3[i] = x3[i] - x2[i];
    }
    
    crossProduct3(c1, b2, b3);
    crossProduct3(c2, b1, b2);
    
    arg1 = dotProduct3(b1, c1);
    b2_norm = sqrt(dotProduct3(b2, b2));
    arg1 = arg1 * b2_norm;

    arg2 = dotProduct3(c2, c1);
    return atan2(arg1, arg2)*180.0 / 3.141592653589793;
}


extern PyObject *
py_dihedral(PyObject *self, PyObject *args)
{
    PyObject *in1, *in2, *in3, *in4;
    PyArrayObject *vec1, *vec2, *vec3, *vec4;
    if (!PyArg_ParseTuple(args, "OOOO", &in1, &in2, &in3, &in4)) return NULL;
    vec1 = (PyArrayObject *) PyArray_ContiguousFromObject(in1, PyArray_DOUBLE, 1, 1);
    if (vec1 == NULL) return NULL;
    vec2 = (PyArrayObject *) PyArray_ContiguousFromObject(in2, PyArray_DOUBLE, 1, 1);
    if (vec2 == NULL) return NULL;
    vec3 = (PyArrayObject *) PyArray_ContiguousFromObject(in3, PyArray_DOUBLE, 1, 1);
    if (vec3 == NULL) return NULL;
    vec4 = (PyArrayObject *) PyArray_ContiguousFromObject(in4, PyArray_DOUBLE, 1, 1);
    if (vec4 == NULL) return NULL;
    
    return PyFloat_FromDouble( dihedral( 
        (double *)vec1->data, (double *)vec2->data,
        (double *)vec3->data, (double *)vec4->data
    ) );
}


/*
 * Bind Python function names to our C functions
 */
static PyMethodDef cgeom_methods[] = {
    {"dihedral", py_dihedral, METH_VARARGS, "calculate dihedral"},
    {NULL, NULL, 0, NULL}
};


/*
 * Python calls this to let us initialize our module
 */
PyMODINIT_FUNC
initcgeom(void)
{
    (void) Py_InitModule("cgeom", cgeom_methods);
    import_array();
}