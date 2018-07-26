//-----------------------------------------------------------------------------------
//
//	Класс трехмерных элементарных частиц
//
//-----------------------------------------------------------------------------------

//#include "tlinear3D.h"
#include "vector3d.h"

//-----------------------------------------------------------------------------------

#ifndef ___Atom3D_h___
#define ___Atom3D_h___

//-----------------------------------------------------------------------------------

class Atom3D // 6*8+8=56 bytes
{
public:
    Atom3D(): v(VECT3D_NULL), n_bonds(0), n_gap(0), color(0)   {}
    Atom3D& operator+=(Vect3d vect)       { v += vect; return *this; }
    Atom3D& operator-=(Vect3d vect)       { v -= vect; return *this; }

    Vect3d R()                     const  { return r; }
    Vector3D& R()                         { return r; }
    Vector3D& V()                         { return v; }
    //int& N_Bonds()                        { return n_bonds; }
    //int& N_Gap()                          { return n_gap; }
    char& N_Bonds()                        { return n_bonds; }
    char& N_Gap()                          { return n_gap; }
    char& Color()                          { return color; }
	double& Init_Z()						   { return init_Z; }

	void Set_Init_Z(double InitZ_new)		{ init_Z = InitZ_new; }
private:
    Vector3D r, v;
    //int n_bonds;                          // количество связей
    //int n_gap;                            // количество разорванных связей
    char n_bonds;                          // количество связей
    char n_gap;                            // количество разорванных связей
    char color;                            // цвет
	double init_Z;
};

//-----------------------------------------------------------------------------------
#include "_list.h"

class Bond3Dv2;
#define Area3d const Area3D&

class Area3D
{
public:
    Area3D() : shearYX(0)                     {}
    Area3D(double x1, double y1, double z1,
           double x2, double y2, double z2)
           : shearYX(0)                       { Set(x1, y1, z1, x2, y2, z2); }
	Area3D(Area3d area)                       { Set(area); }

    void Set(double x1, double y1, double z1,
             double x2, double y2, double z2) { r1.Set(x1, y1, z1); r2.Set(x2, y2, z2); }
    void Set(Area3d area)                     { r1 = area.r1; r2 = area.r2; periodic = area.periodic; shearYX = area.shearYX; }
	void Set(Vect3d r, double ex)             { Set(r, ex, ex, ex); }
	void Set(Vect3d r, double dx,
             double dy, double dz)            { Set(r.x - dx, r.y - dy, r.z - dz, r.x + dx, r.y + dy, r.z + dz); }
    void Set(Vect3d r_1, Vect3d r_2)          { r1 = r_1; r2 = r_2; }

    double X1()                 const   { return r1.x; }
    double Y1()                 const   { return r1.y; }
    double Z1()                 const   { return r1.z; }
    double X2()                 const   { return r2.x; }
    double Y2()                 const   { return r2.y; }
    double Z2()                 const   { return r2.z; }

    Vect3d R1()                 const   { return r1; }
    Vect3d R2()                 const   { return r2; }

    double Width()              const   { return X2() - X1(); }
    double Height()             const   { return Y2() - Y1(); }
    double Depth()              const   { return Z2() - Z1(); }
    double Volume()             const   { return Width()*Height()*Depth(); }

	int IsPeriodic()            const   { return periodic; }
	int&  Periodic()                    { return periodic; }
    double ShearYX()            const   { return shearYX; }

    void ExtendX(Atom3D* as, int n, double scal);
	void ExtendZ(Atom3D* as, int n, double scal);
	void ExtendY(Atom3D* as, int n, double scal);
    void AddShearYX(Atom3D* as, int n, double ShearYX_);  // для Bond3Dv2
	void AddShearYZ(Atom3D* as, int n, double ShearYX_);
	void AddShearZX(Atom3D* as, int n, double ShearYX_);
	void SpinX(Atom3D* as, int n, double spin);
	void SpinY(Atom3D* as, int n, double spin);
	void SpinZ(Atom3D* as, int n, double spin);
	
	

private:
	Vector3D r1, r2;
	int periodic;
    double shearYX;
};

#include "IarrayS.h"

//void SetTypeMaterial(int type_material, double d_min_d, double* pL_pow=0, double* pA_cut=0);
void SetTypeMaterial(double d_min_d, double L_pow=1.0, double a_cut=1.9);
void SetB_dt(double k);
double GetV2(Atom3D* as, int n);
double GetLastAvDisplacement(Atom3D* as, int n);
double GetFullAvDisplacement(Atom3D* as, int n, Vector3D* points0);
double GetFullAvDisplacementGap(Atom3D* as, int n, Vector3D* points0);

void Centrate (Atom3D* as, long n);
void Compaction(Atom3D* as, long n, Area3d area);
double GetDmin_D();
void SetAcut12(double a_cut1, double a_cut2);
double GetAcut();
double GetL_pow();
int GetTypeMaterial();
void CreateMap(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area);
void DetectGroup(Atom3D* as, int n, Array<Bond3Dv2>& bonds, ArrayS<int>& group);



#endif // ___Atom3D_h___


