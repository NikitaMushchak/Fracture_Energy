//-----------------------------------------------------------------------------------
//
//	Saving and loading particle coordinates in A3R format
//
//	Anton M. Krivtsov
//
//	11.04.2001
//  Modified 25.05.2002 by V.Tsaplin
//  Modified 23.05.2011 by V.Tsaplin
//
//-----------------------------------------------------------------------------------

#ifndef ___a3r_h___
#define ___a3r_h___

#include "vector3d.h"

//-----------------------------------------------------------------------------------

struct A3R_HEADER
{
	A3R_HEADER();		// Initialization of the structure members

	char file_type[4];	// "a3r"
	int count;			// number of particles
	int data_start;		// address of the start of the particles data
	char version[10];	// version of a3r format
	double r;			// particle radius
	int count_1;		// reserved for the future use;
};

//-----------------------------------------------------------------------------------

int Save_A3R(const char* file_name, Vector3D* start, const int n, const double r);
Vector3D* Load_A3R(const char* file_name, int& n, double& r);

//-----------------------------------------------------------------------------------

template <class type>
int Save_A3R(const char* file_name, type* start, const int n, const double r)
{
	Vector3D *vect = new Vector3D[n];
	Vector3D *v = vect;

	type* i;

	for (i = start; i < start + n; i++, v++)  { *v = i->R(); }

	Save_A3R(file_name, vect, n, r);

	delete [] vect;

	return 1;
}

//-----------------------------------------------------------------------------------

#include "_list.h"

template <class type>
int Save_A3R(const char* file_name, _List<type>& list, const double r)
{
    int n = list.GetCount();
	Vector3D *vect = new Vector3D[n];
	Vector3D *v = vect;

    int more = list.Restore();

	for (; more; v++)
	{
        type* i = list.Iterate(more);
		*v = i->R();
	}

	Save_A3R(file_name, vect, n, r);

	delete [] vect;

	return 1;
}

#include <string.h>

#endif //___a3r_h___

