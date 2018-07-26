

#include "a3r_b.h"
#include <stdio.h>
#include "colors.h"

//-------------------------------------------------------------------------------------------

A3R_HEADER_B::A3R_HEADER_B() 
: count(0)
, r(0)
, count_1(0)
{ 
	strcpy(file_type, "a3r"); 
	strcpy(version, "b"); 
	data_start = sizeof(A3R_HEADER_B);
    memset(&box, 0, sizeof(float)*6);
}

typedef unsigned char	COLOR_INDEX_C;	// index of the color (short version)
typedef float			COLOR_INDEX;	// index of the color 

//-------------------------------------------------------------------------------------------

int Save_A3R_B(const char* file_name, Atom3D* start, const int n, const double r,
               Box3D_SFLOAT* pBox/*=0*/,  int use_color/*=0*/)
{
	FILE *f=fopen(file_name, "wb");

	A3R_HEADER_B header;
	header.count = n;
	header.r = r;
    if (pBox)  {header.box=*pBox; }
	if (use_color) strcpy(header.version, "d");

  	fwrite(&header, sizeof(header), 1, f);
    Atom3D* i;
    for (i = start; i != start + n; i++) 
    {
        float sr[3] = {(float)i->R().x, (float)i->R().y, (float)i->R().z};
        fwrite(&sr, 3*sizeof(float), 1, f);
    }

    if (use_color)
    {
        for (i = start; i != start + n; i++)		
        {
            COLORREF c = i->Color();
            fwrite(&c, sizeof(COLORREF), 1, f);
        }
    }

	fclose(f);
	return 1;
}

//-----------------------------------------------------------------------------------

int Save_A3R_B(const char* file_name, _List<Atom3D>& list, const double r,
               Box3D_SFLOAT* pBox/*=0*/,  int use_color/*=0*/)
{
	FILE *f=fopen(file_name, "wb");

	A3R_HEADER_B header;
    header.count = list.GetCount();
	header.r = r;
    if (pBox)  {header.box=*pBox; }
	if (use_color) strcpy(header.version, "d");

  	fwrite(&header, sizeof(header), 1, f);
    int more = list.Restore();
	while (more)
	{
        Atom3D* i = list.Iterate(more);
        float sr[3] = {(float)i->R().x, (float)i->R().y, (float)i->R().z};
        fwrite(&sr, 3*sizeof(float), 1, f);
	}

    if (use_color)
    {
        more = list.Restore();
	    while (more)
	    {
            Atom3D* i = list.Iterate(more);
            COLORREF c = i->Color();
            fwrite(&c, sizeof(COLORREF), 1, f);
        }
    }

	fclose(f);
	return 1;
}

//-----------------------------------------------------------------------------------
