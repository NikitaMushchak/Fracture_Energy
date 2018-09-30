// a3r.cpp

#include "a3r.h"
#include "vector3d.h"
#include "util.h"
#include <stdio.h>
//using namespace std;
//-------------------------------------------------------------------------------------------

A3R_HEADER::A3R_HEADER() 
: count(0)
, r(0)
, count_1(0)
{ 
	strcpy(file_type, "a3r"); 
	strcpy(version, "a"); 
	data_start = sizeof(A3R_HEADER);
}

//-------------------------------------------------------------------------------------------

const int nbuf = 40;
const int perm_header[40] = 
{
	0,   1,  2,  3,
	7,   6,  5,  4,
	11, 10,  9,  8,
	12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
	22, 23,
	31, 30, 29, 28, 27, 26, 25, 24,
	35, 34, 33, 32,
	36, 37, 38, 39
};

//-------------------------------------------------------------------------------------------

inline BOOL fopen_m(FILE **file_out, const char *filename, const char *mode)
{
	*file_out = fopen(filename, mode);

	//if ((int)*file_out == 0 || (int)*file_out == -1)
	if (*file_out == 0)
		return 1;
	else
		return 0;
}

int Save_A3R(const char* file_name, Vector3D* start, const int n, const double r) 
{
	FILE *f;
	if (fopen_m(&f, file_name, "wb"))     return 0;


	A3R_HEADER header;
	header.count = n;
	header.r = r;

	if (!is_inverse_byte_order()) // IBM-процессор
	{
		char header_buf[nbuf];
		char *p = reinterpret_cast<char*>(&header);
		for (int k = 0; k < nbuf; k++)
		{
		  header_buf[perm_header[k]] = p[k];
		}
		fwrite(&header_buf, nbuf, 1, f);

		for (Vector3D* i = start; i != start + n; i++) 
		{
            float sr[3] = {(float)i->x, (float)i->y, (float)i->z};
            swap_float(sr[0]);
			swap_float(sr[1]);
            swap_float(sr[2]);

			fwrite(sr, 12, 1, f);
		}
	}
	else
	{  //Intel
		fwrite(&header, sizeof(header), 1, f);
		for (Vector3D* i = start; i != start + n; i++) 
		{
            float sr[3] = {(float)i->x, (float)i->y, (float)i->z};
			fwrite(&sr, 3*sizeof(float), 1, f);
		}
	}
	fclose(f);		

	return 1;
}

//-------------------------------------------------------------------------------------------
// Intel only

Vector3D* Load_A3R(const char* file_name, int& n, double& r)	
{
	FILE *f;
	if (fopen_m(&f, file_name, "rb"))     return 0;

	A3R_HEADER header;

	fread(&header, sizeof(A3R_HEADER),1 , f);
	if(strcmp(header.file_type, "a3r")) return NULL;
	n = header.count;
	r = header.r;

	float* buf = new float[3*n];
	float* j = buf;

	fread(buf, 3*sizeof(float), n, f);
	fclose(f);


	Vector3D* start = new Vector3D[n];
	Vector3D* stop = start + n;
    for (Vector3D* i = start; i != stop; i++)
    { i->x = *j++; i->y = *j++; i->z = *j++; }
	delete [] buf;

	return start;
}
//-------------------------------------------------------------------------------------------



