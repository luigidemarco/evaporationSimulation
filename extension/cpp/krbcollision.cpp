/**
* krbcollsion.cpp
* C++ extension module for finding collision pairs.
**/

#include <Python.h>
#include <vector>
#include <math.h>

typedef struct coords_t {
	double x;
	double y;
} coords_t;

typedef struct pair_coords_t {
	double x;
	double y;
	int i1; // index of particle 1
	int i2; // index of particle 2
} pair_coords_t;

// Find the pairs of particles within distance colcut
// assuming the input vector r is sorted
// by the x-coordinate
static std::vector<pair_coords_t>
find_pairs(std::vector<coords_t> r, double colcut) {
	int l = (int) r.size();
	std::vector<pair_coords_t> out;

	for (int i=0; i<l; i++) {
		for (int j=i+1; j<l; j++) {
			double dx = fabs(r[i].x - r[j].x);
			if (dx < colcut) {
				double dy = fabs(r[i].y - r[j].y);
				if (dy < colcut) {
					pair_coords_t p;
					p.x = dx; p.y = dy;
					p.i1 = i; p.i2 = j;

					out.push_back(p);
				}
			}
			else {
				break;
			}
		}
	}
	return out;
}

// Translates Python objects to C++ objects and vice versa
// C++ extension boilerplate
// See https://docs.python.org/2.7/extending/extending.html
static PyObject*
krbcollision_find_pairs(PyObject* self, PyObject* args) {
	PyObject* input;
	double colcut;

	if (!PyArg_ParseTuple(args, "Od", &input, &colcut))
		return NULL;

	int size = (int) PyList_Size(input);
	std::vector<coords_t> r;
	r.resize(size);

	for (int i=0; i<size; i++) {
		PyObject* o = PyList_GET_ITEM(input, i);
		r[i].x = PyFloat_AS_DOUBLE(PyList_GET_ITEM(o, 0));
		r[i].y = PyFloat_AS_DOUBLE(PyList_GET_ITEM(o, 1));
	}

	std::vector<pair_coords_t> out = find_pairs(r, colcut);
	int l = (int) out.size();

	PyObject* out_pairs = PyList_New(l);
	PyObject* out_coords = PyList_New(l);
	for (int i=0; i<l; i++) {
		pair_coords_t temp = out[i];

		PyObject* op = Py_BuildValue("[ii]", temp.i1, temp.i2);
		PyList_SET_ITEM(out_pairs, i, op);

		PyObject* oc = Py_BuildValue("[dd]", temp.x, temp.y);
		PyList_SET_ITEM(out_coords, i, oc);
	}

	PyObject* output = PyList_New(2);
	PyList_SET_ITEM(output, 0, out_pairs);
	PyList_SET_ITEM(output, 1, out_coords);
	
	return output;
}

// C++ extension boilerplate
// See https://docs.python.org/2.7/extending/extending.html
static PyMethodDef CollisionMethods[] = {
	{"find_pairs", krbcollision_find_pairs, METH_VARARGS,
	"Find collision pair candidates."},
	{NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initkrbcollision(void) {
	(void) Py_InitModule("krbcollision", CollisionMethods);
}
