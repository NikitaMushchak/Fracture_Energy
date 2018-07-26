//-----------------------------------------------------------------------------------
//
//	Saving and loading particle coordinates in A3R format
//
//	Anton M. Krivtsov
//
//	11.04.2001
//	06.10.2002 Modified
//  Modified 19.12.2012 by V.Tsaplin
//  Modified 03.05.2017 by V.Tsaplin
//
//-----------------------------------------------------------------------------------

#ifndef ___a3r_b_h___
#define ___a3r_b_h___

//-----------------------------------------------------------------------------------

#include "atom3D.h"

struct Box3D_SFLOAT     { float r[6]; };

struct A3R_HEADER_B
{
	A3R_HEADER_B();		// Initialization of the structure members

	char file_type[4];	// "a3r"
	long count;			// number of particles
	long data_start;	// address of the start of the particles data
	char version[10];	// version of a3r format
	double r;			// particle radius
	long count_1;		// reserved for the future use;
	Box3D_SFLOAT box;	// box associated with particles;
};

//-------------------------------------------------------------------------------------------

int Save_A3R_B(const char* file_name, Atom3D* start, const int n, const double r,
               Box3D_SFLOAT* pBox = 0,  int use_color = 0);

int Save_A3R_B(const char* file_name, _List<Atom3D>& list, const double r,
               Box3D_SFLOAT* pBox = 0,  int use_color = 0);

//-----------------------------------------------------------------------------------

//Atom3D* Load_A3R(const char* file_name, long& n, double& r, Box3D* box = NULL);	

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

#include <string.h>

#endif //___a3r_b_h___
