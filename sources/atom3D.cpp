#include "atom3D.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdio.h>

///////////////////////////////////////////////////////////////////////////////////////

Area3D GetBoundBox(Atom3D* as, long n)
{
    Vector3D r1 = as->R();
    Vector3D r2 = as->R();
    Atom3D* i = as+1;
	Atom3D* i_end = as + n;
    for (; i != i_end; i++)
    {
        if (i->R().x < r1.x) { r1.x = i->R().x; }
        if (i->R().y < r1.y) { r1.y = i->R().y; }
        if (i->R().z < r1.z) { r1.z = i->R().z; }
        if (i->R().x > r2.x) { r2.x = i->R().x; }
        if (i->R().y > r2.y) { r2.y = i->R().y; }
        if (i->R().z > r2.z) { r2.z = i->R().z; }
    }
    Area3D area; area.Set(r1, r2);
    return area;
}

void Centrate (Atom3D* as, long n)
{
  /*  Vector3D rc = VECT3D_NULL;
    Atom3D* i = as;
	Atom3D* i_end = as + n;
    for (; i != i_end; i++) rc += i->R();
    rc /= n;
    for (i = as; i != i_end; i++) i->R() -= rc; */

    Area3D area = GetBoundBox(as, n);
    Vector3D r1 = area.R1();
    Vector3D r2 = area.R2();
    r1+=r2; r1/=2;
    Atom3D* i = as;
	Atom3D* i_end = as + n;
    for (; i != i_end; i++)   { i->R() -= r1; }
}

void Compaction(Atom3D* as, long n, Area3d area)
{
    double	  x1 = area.X1(),     y1 = area.Y1(),     z1 = area.Z1();
    double	  x2 = area.X2(),     y2 = area.Y2(),     z2 = area.Z2();
	double    ax = area.Width(),  ay = area.Height(), az = area.Depth();

    int m;
    for (m = 0; m < n; m++)
    {
        if      (as[m].R().x >= x2)   as[m].R().x -= ax;
        else if (as[m].R().x <  x1)   as[m].R().x += ax;
        if      (as[m].R().y >= y2)   as[m].R().y -= ay;
        else if (as[m].R().y <  y1)   as[m].R().y += ay;
        if      (as[m].R().z >= z2)   as[m].R().z -= az;
        else if (as[m].R().z <  z1)   as[m].R().z += az;
    }
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

double GetV2(Atom3D* as, int n)
{
    double K = 0.0;
    Atom3D* i;
    Atom3D* i_end = as + n;

    for (i = as; i != i_end; i++)
    {
        K += i->V().Sqr();
    }

    return K;
}

/////////////////////////////////////////////////////////////////////////////////
//static int SecBonds=0;
//static int nSecBonds=0;

//////////////////////////////////////////////////////////////////////////////


