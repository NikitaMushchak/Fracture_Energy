//-----------------------------------------------------------------------------------
//
//	Класс трехмерных элементарных частиц, версия 2
//
//-----------------------------------------------------------------------------------

#include "atom3D.h"
#include "tlinear3D.h"

//-----------------------------------------------------------------------------------

#ifndef ___Atom3Dv2_h___
#define ___Atom3Dv2_h___

//-----------------------------------------------------------------------------------

// связи с возможностью добавления угловых пружинок
class Bond3Dv2 // 8+8+3*8+8=48
{
public:
    Bond3Dv2(): shift(0), gap(1)                               {}
    Bond3Dv2(int i0, int j0, char sh=0, float eps_L = 0.0003):
           i(i0), j(j0), shift(sh), gap(1), eps_L(0.0003)                     {}

    void FixLength(Atom3D* as, Area3d area, double d_av);
    void FixLength(Atom3D* as, Area3d area, double d_av, double C_bond);
    int First()              const { return i; }        // индексы массива частиц Atom3D
    int Second()             const { return j; }
    int Shifted()            const { return shift; }
    Vector3D& dR()                 { return _dR; }
    double Length()          const { return length; }
	double InitZ1()          const { return initZ1; }
	double InitZ2()          const { return initZ2; }
    double C_bond()          const { return C; }
	double Eps_L()          const { return eps_L; }
    Vector3D& dR(Atom3D* as, Area3d area);              // dR вычисляется заново
    float& ForceL()                { return force_l; }
    char& Gap()                    { return gap; }
    static char& GapM();
	void Set_C(double C_new)		{ C = C_new; }
	void Set_Eps(float Eps_new) { eps_L = Eps_new; }
	void Set_InitZ1(float InitZ1_new) { initZ1 = InitZ1_new; }
	void Set_InitZ2(float InitZ2_new) { initZ2 = InitZ2_new; }

private:
    int i;
	int j;
    char shift;        // связь между частицами соседних периодических областей
    char gap;          // ==1  - разрываемая связь, ==2  - разорванная связь
    float length;      // длина связи при нулевой силе взаимодействия
    Vector3D _dR;      // вектор R2-R1 с учетом смещения периодических областей
    float force_l;     // сила взаимодействия
    float C;           // коэф. упругости (между F и delta L)
	float eps_L;		//максимальное допустимое удлинение связи
	float initZ1;
	float initZ2;
};

// угловые пружинки
class AngleBond3D // 4*5=20
{
public:
    AngleBond3D()                                              {}
    AngleBond3D(int i0, int j0): i(i0), j(j0)                  {}
    AngleBond3D(int i0, int j0, int k0): i(i0), j(j0), k(k0)   {}

    void FixAngle(Bond3Dv2* bonds, double d_av);
    void FixAngle(Bond3Dv2* bonds, double d_av, double C_bond);
    int First()              const { return i; } // индексы массива связей Bond3Dv2
    int Second()             const { return j; } // (списки Bond3Dv2 использоваться не будут)
    int Third()              const { return k; }
    double Angle()           const { return angle; }
    double GetC()            const { return C; }
	void SetC(double C_new) { C = C_new; }
private:
    int i, j, k;
    float angle;      // угол между связями при нулевой деформации
    float C;          // угловая жесткость
};

void CreateConfig3D(Array<Atom3D>& particles, int n_layers, Area3D& area, double d_av,
                    Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, double slope=0);

void SaveConfig(char* file_name, Array<Atom3D>& particles, Area3D& area,
                Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, int& s, double* V2_max=0);
void LoadConfig(char* file_name, Array<Atom3D>& particles, Area3D& area,
                Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, int& s, double* V2_max=0);
void SaveConfig2(char* file_name, Array<Vector3D>& points, int N_equ);
void LoadConfig2(char* file_name, Array<Vector3D>& points, int& N_equ);

void Step3D(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds, int p,
	double press = 0.0, Vector3D normal = { 0.,0.,0. }, Vector3D center = {0.,0.,0.},double diameter = 0. , double* pFacet=0, _List<Atom3D>* pListSlidingSealingY=0);

void GetStress_XYZ(Area3d area, Atom3D* as, Array<Bond3Dv2>& bonds, double &(sigma_x), double &(sigma_y), double &(sigma_z), double &(sigma_xy), double &(sigma_xz), double &(sigma_yz));

double GetNBonds(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, int axis, double section);

Tensor GetStressTensor(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, double &(sigma_x), double &(sigma_y), double &(sigma_z), double &(sigma_xy), double &(sigma_xz), double &(sigma_yz));
Vector3D GetStress(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m,
                   int axis, double section); // axis = 1,2,3;  0<section<1

void SetAngleSprings(double C_fi, double a_cut1=0, double a_cut2=0);
double GetC_fi();

void CreateFracture(Atom3D* as, Area3d area, Array<Bond3Dv2>& bonds,
                    Vector3D center, double diameter, Vector3D normal);

void Strain_XYZ_Gap(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds);
void Stress_XYZ_Gap(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds, int num);
void Get_True_Strain(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds);
void Get_Init_Coord(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds);
void Get_Press(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int n, int m, AngleBond3D* a_bonds, double d_S, double &(P));
void Get_True_Strain_N(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds, double P, double R);
void Get_Init_Coord_N(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds);

int Save_A3R_Gap(const char* file_name, Atom3D* as, int n, Array<Bond3Dv2>& bonds, Area3d area, const double r);


//double Volume(Atom3D* as, Area3d area, Array<Bond3Dv2>& bonds);
double Volume(Array<Bond3Dv2>& bonds,double ds=1);
double Energy_Potencial(Array<Bond3Dv2>& bonds,double ds=1);
int GetLiq(Array<Bond3Dv2>& bonds);

double GetFullAvDisplacement(Atom3D* as, int n, Atom3D* points0);

void CreateConfig3D(Array<Atom3D>& particles,
                    int n_layersX,
                    int n_layersZ,
                    int n_layersY1,
                    int n_layersY2,
                    int n_layersY3,
                    Area3D& area, double d_av,
                    Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds);
/*Aspect ratio для актуального состояния трещины*/
double AspectRatio(Atom3D* as, Array<Bond3Dv2>& bonds);

/*Создание различных упругих модулей*/
/*factor - множитель, на который изменяются жесткости связей, part - часть для которой изменяется жесткость (-1 левая часть, 0 - правая часть, 1 - правая часть, все в отношении вдоль оси Х)*/
void SetDifferentModules(Area3D& area, Atom3D* as, Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, double factor, int part,double part_third_centr);

/*Создание различных трещиностойкостей*/
/*factor - множитель, на который изменяются критическое удлинение связей, part - часть для которой изменяется жесткость (-1 левая часть, 0 - правая часть, 1 - правая часть, все в отношении вдоль оси Х)*/
void SetDifferentToughness(Area3D& area, Atom3D* as, Array<Bond3Dv2>& bonds, double factor, int part, double part_third_centr);
void Save_XYZ_Gap(Area3D & area, int n_step, Atom3D* as, Array<Bond3Dv2>& bonds, int n,double dV);
void Save_XYZ_Media(Area3D & area, int n_step, Atom3D * as, int n, double E_centr, double E_left, double E_right, double KC_centr, double KC_left, double KC_right, double part_third);



#endif // ___Atom3Dv2_h___


