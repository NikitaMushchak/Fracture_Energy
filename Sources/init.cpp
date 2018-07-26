#include "atom3Dv2.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <cmath>

static const double C_bond = 1;
static double L_pow = 1;

/////////////////////////////////////////////////////////////////////////////////

static int type_material=2;//0,1,2,3
static double d_min_d=1;
static double a_cut = 1.9;

double GetDmin_D()    { return d_min_d; }
int GetTypeMaterial() { return type_material; }
double GetAcut()      { return a_cut; }
double GetL_pow()     { return L_pow; }

void Set_d_min_d_step(double d_min_d, double L_pow);

void SetTypeMaterial(double d_min_d, double L_pow/*=1*/, double A_cut/*=1.9*/)
{
    ::d_min_d=d_min_d;
    ::L_pow = L_pow;
    a_cut = A_cut;

    Set_d_min_d_step(d_min_d, L_pow);
}

static void _12romb_cube(Vector3D& r1)
{
    // область r1: заполняется часть ромбододекаэдра, исключая вписанный в него куб
    r1.SetRand(-0.25,-0.25,-0.25,0.25,0.25,0.25);
    double mxy = max2(fabs(r1.x),fabs(r1.y));
	if (fabs(r1.z) >= mxy) {
		r1.z -= 0.5* abs(r1.z) / r1.z;;
	}//_copysign(0.5, r1.z); }
	else if (fabs(r1.x) > fabs(r1.y)) { 
		r1.x -= 0.5* abs(r1.x) / r1.x;;
	}//_copysign(0.5, r1.x); }
	else {
			r1.y -= 0.5* abs(r1.y) / r1.y;
	}//_copysign(0.5, r1.y); } }
}

int Ceil(int a, int b) // b > 0
{
    int A = 0;
    if (a > 0) { A = a / b + 1; }
    return A - (A*b - a) / b;
}


// создание материала на основе ГЦК
void NewSystem3D_02(Array<Atom3D>& particles, int n_layers, Area3D& area, double d_av, double slope)
// угол наклона 0<=slope<PI/2
{
    // для изотропии type_material=1 -> d_min=0.22*d_av; (при a_cut=1.3*d_av)
    // для изотропии type_material=2 -> d_min=0.30*d_av;
    double eps = d_min_d;
    double dp = d_av * M_SQRT2;         //  толщина слоя частиц (двойного)
    double dp_2 = dp / 2;
    //double D =   dp   * n_layers;       //  размер области (по оси Z)
    double D_2 = dp_2 * n_layers;
    int N2 = n_layers * 2;              //  количество одинарных слоев (по оси Z)

    // n11+n12 должно быть четным
    double Nc_fi=n_layers*cos(slope);
    double Ns_fi=n_layers*sin(slope);
    int m1=(int)floor(Nc_fi+Ns_fi+0.5);
    int m2=(int)floor(Nc_fi-Ns_fi+0.5);
    int n11=m1+m2;
    int n12=m1-m2;
    int n2=n11*n11+n12*n12;
    int n = n2*n_layers;                //  количество частиц

    double D_XY = dp_2 * sqrt(float(n2));      //  размер области (по осям X,Y)
    double D_2XY = D_XY / 2;
    area.Set(VECT3D_NULL, D_2XY, D_2XY, D_2);
    particles.Create(n);

    //double z1 = dp/4 - D_2;
    Vector3D turn;    turn.SetTurn(-atan2(float(n12), float(n11)), Vector3D(0,0,1));

    int i, j, k, m = 0;
    for(k = 0;  k < N2; k++)
    for(j = 0;	j < n11+n12; j++)
    {
        // диапазон индексов i находится согласно границам повернутого квадрата на угол slope
        //int i_min = (int)ceil(j<=n11 ? -j*n12/(double)n11 : (j*n11-n2)/(double)n12);
        //int i_max = (int)ceil(j<n12  ?  j*n11/(double)n12 : (n2-j*n12)/(double)n11);
        int i_min = j<=n11 ? Ceil(-j*n12, n11) : Ceil(j*n11-n2, n12);
        int i_max = j<n12  ? Ceil( j*n11, n12) : Ceil(n2-j*n12, n11);
        if ((i_min^j^k)&1)   { i_min++; }      // четность индексов узлов ГЦК
        for(i = i_min;	i < i_max; i+=2)
        {
            Vector3D r1;

            if (type_material==0)                                           // ромбододекаэдр
            {
                r1.SetRand(-0.25,-0.25,-0.25,0.25,0.25,0.75);
                if (r1.z > 0.25)  { _12romb_cube(r1); }
            }
            else { if (type_material==1)  { _12romb_cube(r1); }             // ромбодод. - куб
            else { do  { _12romb_cube(r1); }  while (r1.Sqr()<0.125); } }   // ромбодод. - куб - сфера

            r1 *= dp * (1 - eps);
            r1.x += i*dp_2;
            r1.y += j*dp_2;
            r1.z += k*dp_2;
            r1.Rotate(turn, r1);
            //r1.x -= D_2XY;//r1.y -= D_2XY;//r1.z += z1;

            particles[m++].R() = r1;
        }
    }

    if (m!=n) // ошибка
    {
        printf("m=%i\tn=%i\n", m, n);
        exit(0);
    }

    Atom3D* as = particles.GetData();
    Centrate(as, n);
    Compaction(as, n, area);
}

//static const double M_SQRT3=sqrt(3);

/*
// случайное накидывание частиц
void NewSystem3D_5(Array<Atom3D>& particles, int n_layers, Area3D& area, double d_av)
{
}*/

void NewSystem3D(Array<Atom3D>& particles, int n_layers, Area3D& area, double d_av, double slope/*=0*/)
{
    //if (type_material==3)           { NewSystem3D_3 (particles, n_layers, area, d_av, slope); }

    if (type_material==5)   { /*NewSystem3D_5 (particles, n_layers, area, d_av);*/ }
    else                    { NewSystem3D_02(particles, n_layers, area, d_av, slope); }
}

// сортировка частиц по ячейкам для обнаружения соседних частиц
void CreateMap(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area)
{
	double		ax = area.Width(),		ay = area.Height(),		az = area.Depth();
	double		x0 = area.X1(),			y0 = area.Y1(),			z0 = area.Z1();

	// --- Вычисление пространственных параметров ---

	Atom3D* i;
	Atom3D* i_end = as + n;
	long ix, iy;
	int iz;

	long	nx = long(ax / a_cut),	ny = long(ay / a_cut),	nz = long(az / a_cut);

	nx = nx < 1 ? 1 : nx; ny = ny < 1 ? 1 : ny; nz = nz < 1 ? 1 : nz;

	double		cx = ax / nx,			cy = ay / ny,			cz = az / nz;

	cells.Create(nx, ny, nz);

	for (i = as; i != i_end; i++)
	{
		ix = long((i->R().x - x0) / cx);	if(ix < 0 || ix >= nx) continue;
		iy = long((i->R().y - y0) / cy);	if(iy < 0 || iy >= ny) continue;
		iz = long((i->R().z - z0) / cz);	if(iz < 0 || iz >= nz) continue;

		cells[ix][iy][iz].Insert(i);
	}
}

//////////////////////////////////////////////////////////////////////////////

static char gap_m=2;//50;//2
char& Bond3Dv2::GapM()  { return gap_m; }

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//  связи с угловыми пружинками

Vector3D& Bond3Dv2::dR(Atom3D* as, Area3d area)
{
    _dR.Difference(as[j].R(), as[i].R());

    // изменение радиус-вектора второй частицы относительно первой
    // с учетом того, что эти частицы находятся в соседних периодических областях

    if (shift)
    {
        // если частицы i, j находятся в разных периодических областях,
        // вычисляем соответствующие поправки Shift к вектору dR
        Vector3D Shift(0,0,0);

        //if (shift&3)  { Shift.x += (shift&1)  ? -area.Width()  : area.Width();  }
        // не бывает ((shift&2) != 0), см. DetectBonds1(...,_List<Bond3Dv2>&)
        if (shift&1)  { Shift.x -= area.Width();  }
        if (shift&12) { Shift.y += (shift&4)  ? -area.Height() : area.Height(); }
        if (shift&48) { Shift.z += (shift&16) ? -area.Depth()  : area.Depth();  }

        // поправка вектора периода с учетом деформации сдвига
        // периодической области

        double shearYX = area.ShearYX();
        if (shearYX)  { Shift.x += Shift.y * shearYX; }

        _dR += Shift;
    }
    return _dR;
}

//  Фиксируем длину связи (принимаем текущую длину связи за длину недеформированной связи)
void Bond3Dv2::FixLength(Atom3D* as, Area3d area, double d_av)
{
    length = float(dR(as, area).Abs());             //  длина недеформированной связи
    C=float(::C_bond*pow(d_av/length, L_pow));      //  обновляем жесткость связи
}

void AngleBond3D::FixAngle(Bond3Dv2* bonds, double d_av)
{
    Bond3Dv2& b1 = bonds[First()];
    Bond3Dv2& b2 = bonds[Second()];
    Vector3D& dR1 = b1.dR();
    Vector3D& dR2 = b2.dR();
    double cos_fi =
        (b2.First()==b1.First() || b2.Second()==b1.Second()) ?
        (dR1*dR2) / dR1.Abs() / dR2.Abs():
       -(dR1*dR2) / dR1.Abs() / dR2.Abs();

    angle = float(acos(cos_fi));
    C=float(GetC_fi()*min2(dR1.Sqr(), dR2.Sqr()) / d_av / d_av);  // alpha=2
    //C=float(GetC_fi()*min2(dR1.Abs(), dR2.Abs()) / d_av);       // alpha=1

    /*Vector3D dR3;
    if (b2.First()==b1.First() || b2.Second()==b1.Second())
    {
        dR3.Difference(dR1, dR2);
    }
    else
    {
        dR3.Sum(dR1, dR2);
    }

    C=GetC_fi()*min3(dR1.Sqr(), dR2.Sqr(), dR3.Sqr());*/
}

void Area3D::AddShearYX(Atom3D* as, int n, double ShearYX_)
{
	shearYX = ShearYX_;
	Atom3D* i;
	Atom3D* i_end = as + n;

	// сдвигаем частицы
	for (i = as; i != i_end; i++) { i->R().x += i->R().y * shearYX; }
	//for (i = as; i != i_end; i++) { i->R().y += i->R().x * shearYX; }
}

void Area3D::AddShearZX(Atom3D* as, int n, double ShearYX_)
{
	shearYX = ShearYX_;
	Atom3D* i;
	Atom3D* i_end = as + n;

	// сдвигаем частицы
	for (i = as; i != i_end; i++) { i->R().x += i->R().z * shearYX; }
	for (i = as; i != i_end; i++) { i->R().z += i->R().x * shearYX; }
}

void Area3D::AddShearYZ(Atom3D* as, int n, double ShearYX_)
{
	shearYX = ShearYX_;
	Atom3D* i;
	Atom3D* i_end = as + n;

	// сдвигаем частицы
	for (i = as; i != i_end; i++) { i->R().z += i->R().y * shearYX; }
	for (i = as; i != i_end; i++) { i->R().y += i->R().z * shearYX; }
}

void Area3D::ExtendX(Atom3D* as, int n, double scal)
{
	r1.x *= scal;
	r2.x *= scal;
	Atom3D* i;
	Atom3D* i_end = as + n;

	for (i = as; i != i_end; i++) { i->R().x *= scal; }
}

void Area3D::ExtendZ(Atom3D* as, int n, double scal)
{
	r1.z *= scal;
	r2.z *= scal;
	Atom3D* i;
	Atom3D* i_end = as + n;

	for (i = as; i != i_end; i++) { i->R().z *= scal; }
}

void Area3D::ExtendY(Atom3D* as, int n, double scal)
{
	r1.y *= scal;
	r2.y *= scal;
	Atom3D* i;
	Atom3D* i_end = as + n;

	for (i = as; i != i_end; i++) { i->R().y *= scal; }
}

void Area3D::SpinX(Atom3D* as, int n, double spin)
{
	Atom3D* i;
	Atom3D* i_end = as + n;

	for (i = as; i != i_end; i++) { i->R().y = i->R().y * cos(spin) - i->R().z * sin(spin); }
	for (i = as; i != i_end; i++) { i->R().z = i->R().y * sin(spin) + i->R().z * cos(spin); }
}

void Area3D::SpinY(Atom3D* as, int n, double spin)
{
	Atom3D* i;
	Atom3D* i_end = as + n;

	for (i = as; i != i_end; i++) { i->R().z = i->R().z * cos(spin) - i->R().x * sin(spin); }
	for (i = as; i != i_end; i++) { i->R().x = i->R().z * sin(spin) + i->R().x * cos(spin); }
}

void Area3D::SpinZ(Atom3D* as, int n, double spin)
{
	Atom3D* i;
	Atom3D* i_end = as + n;

	for (i = as; i != i_end; i++) { i->R().x = i->R().x * cos(spin) - i->R().y * sin(spin); }
	for (i = as; i != i_end; i++) { i->R().y = i->R().x * sin(spin) + i->R().y * cos(spin); }
}
//////////////////////////////////////////////////////////////////////////////
// нахождение связей между частицами (связи с угловыми пружинками),
// окончательная версия (используются треугольники из связей)

// as[n]       - массив частиц
// cells[][][] - трехмерный массив ячеек со списками частиц в каждой ячейке
// area        - область 3D простраства
// bonds[]     - массив связей между частицами
// group[i][]  - массив связей при каждой частице i

void DetectBonds0(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, IArray<Bond3Dv2>& bonds, int& m);
void DetectBonds1(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, IArray<Bond3Dv2>& bonds, int& m);

void DetectBonds(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, Array<Bond3Dv2>& bonds, double d_av)
{
    IArray<Bond3Dv2> ibonds;    //  используем бесконечный массив
    int m;
    if (area.IsPeriodic())    DetectBonds1(as, n, cells, area, ibonds, m);
    else                      DetectBonds0(as, n, cells, area, ibonds, m);

    bonds.Create(m);
    int i;
    for (i=0; i<m; i++)
    {
        bonds[i] = ibonds[i];
        bonds[i].FixLength(as, area, d_av);
    }
}

// создание двумерного массива group для
// определения индекса связи по индексу частицы и локальному номеру связи,
// group[индекс частицы][номер связи при этой частице] == индекс этой связи в массиве bonds

void DetectGroup(Atom3D* as, int n, Array<Bond3Dv2>& bonds, ArrayS<int>& group)
{
    group.Create(n, 0);
    long i, j;
    for (i=0; i<n; i++)  { group[i].Create(as[i].N_Bonds()); }

    Array<int> count_group(n); // номер связи при каждой частице
    count_group = 0;

    int m = bonds.GetCount();
    for (i=0; i<m; i++)
    {
        j = bonds[i].First();      group[j][count_group[j]++] = i;
        j = bonds[i].Second();     group[j][count_group[j]++] = i;
    }
}

// связи в непериодической области:
void DetectBonds0(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, IArray<Bond3Dv2>& bonds, int& m)
{
    double		ax = area.Width(),		ay = area.Height(),		az = area.Depth();
    double		x0 = area.X1(),			y0 = area.Y1(),			z0 = area.Z1();

    Atom3D* i;
    Atom3D* i_end = as + n;
    long ix, iy;
    int iz;

   	double aa_cut = a_cut * a_cut;

    int nx = cells.GetCount(),   ny = cells[0L].GetCount(),   nz = cells[0L][0L].GetCount();
	double       cx = ax / nx,   cy = ay / ny,                cz = az / nz;

    m = 0;
    for (i = as; i != i_end; i++)
    {
        ix = long((i->R().x - x0) / cx);
        iy = long((i->R().y - y0) / cy);
        iz = long((i->R().z - z0) / cz);

        // ищем ближайшие частицы в этой и соседних ячейках cells
        long jx, jy;
        int jz;
        int _ix = ix - 1, _iy = iy - 1, _iz = iz - 1;
        int iy_ = iy + 2, iz_ = iz + 2;

        // на границе области свободные края:
        // выходящие за пределы массива cells индексы не рассматриваем
        _ix = _ix < 0 ? 0 : _ix;
        _iy = _iy < 0 ? 0 : _iy;     _iz = _iz < 0 ? 0 : _iz;
        iy_ = iy_ >= ny ? ny : iy_;  iz_ = iz_ >= nz ? nz : iz_;

		for(jx = _ix;	jx <= ix;	jx++) // не рассматриваем соседние ячейки cells справа (исключаем дублирование связей - частиц i,j и j,i)
		for(jy = _iy;	jy < iy_;	jy++)
		for(jz = _iz;	jz < iz_;	jz++)
		{
			_List<Atom3D>* lj = &cells[jx][jy][jz];
			BOOL more = lj->Restore();
			while(more)
			{
				Atom3D* j = lj->Iterate(more);

                if (i == j) { jy = iy_; jz = iz_; break; }  // выходим из трех циклов (по jx, jy, jz) -
                                                            // исключаем дублирование связей - частиц i,j и j,i
				Vector3D dr;
				dr.Difference(j->R(), i->R());
                double dd = dr.Sqr();

                if (dd < aa_cut)
                {
                    bonds[m++] = Bond3Dv2(int(i-as), int(j-as));
                    i->N_Bonds()++;
                    j->N_Bonds()++;
                }
            }
        }
    }
}

// учитываются связи между частицами соседних периодов (условие периодичности):
void DetectBonds1(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, IArray<Bond3Dv2>& bonds, int& m)
{
    double		ax = area.Width(),		ay = area.Height(),		az = area.Depth();
    double		x0 = area.X1(),			y0 = area.Y1(),			z0 = area.Z1();

    Atom3D* i;
    Atom3D* i_end = as + n;
    long ix, iy;
    int iz;

   	double aa_cut = a_cut * a_cut;

    int nx = cells.GetCount(),   ny = cells[0L].GetCount(),   nz = cells[0L][0L].GetCount();
	double       cx = ax / nx,   cy = ay / ny,                cz = az / nz;

    m = 0;
    for (i = as; i != i_end; i++)
    {
        ix = long((i->R().x - x0) / cx);
        iy = long((i->R().y - y0) / cy);
        iz = long((i->R().z - z0) / cz);

        // ищем ближайшие частицы в этой и соседних ячейках cells
        int _jx, _jy, _jz;
        int _ix = ix - 1, _iy = iy - 1, _iz = iz - 1;
        int iy_ = iy + 2, iz_ = iz + 2;

		for(_jx = _ix;	_jx <= ix;	_jx++) // не рассматриваем соседние ячейки cells справа (исключаем дублирование связей - частиц i,j и j,i)
		for(_jy = _iy;	_jy < iy_;	_jy++)
		for(_jz = _iz;	_jz < iz_;	_jz++)
		{
            char shift = 0;
            long jx = _jx, jy = _jy;
            int jz = _jz;

            // ставим периодические граничные условия по всей границе области
            if (jx < 0)    { jx = nx - 1;  shift = 1; }
            // else if (jx == nx)  { jx = 0;  shift = 2; }   (вообще jx <= ix < nx)
            if (jy < 0)    { jy = ny - 1;  shift |= 4; }
            else if (jy == ny)  { jy = 0;  shift |= 8; }
            if (jz < 0)    { jz = nz - 1;  shift |= 16; }
            else if (jz == nz)  { jz = 0;  shift |= 32; }

			_List<Atom3D>* lj = &cells[jx][jy][jz];
			BOOL more = lj->Restore();
			while(more)
			{
				Atom3D* j = lj->Iterate(more);

                if (i == j) { _jy = iy_; _jz = iz_; break; }    // выходим из трех циклов (по _jx, _jy, _jz) -
                                                                // исключаем дублирование связей - частиц i,j и j,i
                Bond3Dv2 bond(int(i-as), int(j-as), shift);
                double dd = bond.dR(as, area).Sqr();

                if (dd < aa_cut)
                {
                    bonds[m++] = bond;
                    i->N_Bonds()++;
                    j->N_Bonds()++;
                }
            }
        }
    }
}

double a_cut1;
double a_cut2;

// диапазон длин связей для установки угловых пружинок (a_cut2 <= a_cut1)
void SetAcut12(double a_cut1, double a_cut2)
{
    ::a_cut1 = a_cut1 ? a_cut1 : a_cut;
    ::a_cut2 = a_cut2 ? a_cut2 : 1.5;
}

// нахождение углов между связями для установки угловых пружинок:
void DetectAngleBonds1(int n, Bond3Dv2* bonds, ArrayS<int>& group,
                      Array<AngleBond3D>& a_bonds)
{
    IArray<AngleBond3D> ia_bonds;    //  используем бесконечный массив
    long i;
    int j, k, t, na_bonds=0;

    for (i = 0; i < n; i++) // все частицы
    {
        int m = group[i].GetCount();

        for (j = 0; j < m; j++)     // все пары связей ...
        {
            Bond3Dv2 &b1 = bonds[group[i][j]];
            if (b1.Length() >= a_cut1 || b1.Length() < a_cut2)    continue;
            long p = (b1.First() == i) ? b1.Second() : b1.First();

        for (k = 0; k < j; k++)     // ... при этой частице
        {
            Bond3Dv2 &b2 = bonds[group[i][k]];
            if (b2.Length() >= a_cut1 || b2.Length() < a_cut2)    continue;
            int  s = (b2.First() == i) ? b2.Second() : b2.First();

            for (t = 0; t < group[p].GetCount(); t++) // все связи при частице as[p]
            {
                Bond3Dv2 &b3 = bonds[group[p][t]];
                int u = (b3.First() == p) ? b3.Second() : b3.First();

                if (u == s)  // между частицами as[p] и as[s] есть связь,
                {            // тогда устанавливаем угловую пружинку
                    if (b3.Length() < a_cut1 && b3.Length() >= a_cut2) {
                        ia_bonds[na_bonds++] = AngleBond3D(group[i][j], group[i][k], group[p][t]); }
                    break;
                }
            }
        }}
    }

    if (na_bonds==0)  return;

    a_bonds.Create(na_bonds);
    for (j = 0; j < na_bonds; j++)  { a_bonds[j] = ia_bonds[j]; }
}

void DetectAngleBonds(int n, Bond3Dv2* bonds, ArrayS<int>& group,
                      Array<AngleBond3D>& a_bonds, double d_av)
{
    DetectAngleBonds1(n, bonds, group, a_bonds);

    int j, na_bonds = a_bonds.GetCount();
    for (j = 0; j < na_bonds; j++)  { a_bonds[j].FixAngle(bonds, d_av);}
}

//////////////////////////////////////////////////////////////////////////////

void CreateConfig3D(Array<Atom3D>& particles, int n_layers, Area3D& area, double d_av,
                    Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, double slope/*=0*/)
{
    //  массив частиц со случайными отклонениями от плотноупакованной решетки
    NewSystem3D(particles, n_layers, area, d_av, slope);

    Atom3D* as = particles.GetData();
    int n = particles.GetCount();
    ArrayS< _List<Atom3D> > cells;

    //  сортировка частиц по ячейкам
    CreateMap(as, n, cells, area);

    //  нахождение массива связей
    DetectBonds(as, n, cells, area, bonds, d_av);

    if (GetC_fi()!=0.0) // если угловая жесткость != 0, устанавливаем угловые пружинки
    {
        //  находим группы связей при каждой частице
        ArrayS<int> group;
        DetectGroup(as, n, bonds, group);

        //  нахождение массива углов между связями
        DetectAngleBonds(n, bonds.GetData(), group, a_bonds, d_av);
    }
}

//////////////////////////////////////////////////////////////////////////////

