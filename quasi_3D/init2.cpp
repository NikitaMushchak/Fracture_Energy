#include "atom3Dv2.h"
#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdio.h>
//
//static void _12romb_cube(Vector3D& r1)
//{
//    // область r1: заполн€етс€ часть ромбододекаэдра, исключа€ вписанный в него куб
//    r1.SetRand(-0.25,-0.25,-0.25,0.25,0.25,0.25);
//    double mxy = max2(fabs(r1.x),fabs(r1.y));
//    if (fabs(r1.z) >= mxy)              { r1.z -= _copysign(0.5, r1.z); }
//    else { if (fabs(r1.x) > fabs(r1.y)) { r1.x -= _copysign(0.5, r1.x); }
//    else                                { r1.y -= _copysign(0.5, r1.y); } }
//}
//
//// материал на основе √÷  (без наклона)
//void NewSystem3D_01(Array<Atom3D>& particles,
//                    int n_layersX, int n_layersY, int n_layersZ, Area3D& area, double d_av)
//{
//    // дл€ изотропии type_material=1 -> d_min=0.22*d_av; (при a_cut=1.3*d_av)
//    // дл€ изотропии type_material=2 -> d_min=0.30*d_av;
//    double eps = GetDmin_D();
//    double dp = d_av * M_SQRT2;         //  толщина сло€ частиц
//    double dp_2 = dp / 2;
//    //double D =   dp   * n_layers;       //  размер области
//    double DX_2 = dp_2 * n_layersX;
//    double DY_2 = dp_2 * n_layersY;
//    double DZ_2 = dp_2 * n_layersZ;
//    int NX2 = n_layersX * 2;
//    int NY2 = n_layersY * 2;
//    int NZ2 = n_layersZ * 2;
//    int n = NX2*NY2*n_layersZ;             //  количество частиц
//
//    area.Set(VECT3D_NULL, DX_2, DY_2, DZ_2);
//    particles.Create(n);
//
//    double x1 = dp/4 - DX_2;
//    double y1 = dp/4 - DY_2;
//    double z1 = dp/4 - DZ_2;
//
//    int type_material=GetTypeMaterial();
//
//    int i, j, k, m = 0;
//    for(i = 0;	i < NX2; i++)
//    for(j = 0;	j < NY2; j++)
//    for(k = (i^j)&1;  k < NZ2; k+=2)
//    {
//        Vector3D r1;
//
//        if (type_material==0)
//        {
//            r1.SetRand(-0.25,-0.25,-0.25,0.25,0.25,0.75);
//            if (r1.z > 0.25)  { _12romb_cube(r1); }
//        }
//        else { if (type_material==1)  { _12romb_cube(r1); }
//        else { do  { _12romb_cube(r1); }  while (r1.Sqr()<0.125); } }
//
//        r1 *= dp * (1 - eps);
//        r1.x += i*dp_2 + x1;
//        r1.y += j*dp_2 + y1;
//        r1.z += k*dp_2 + z1;
//
//        particles[m++].R() = r1;
//    }
//
//    if (m!=n) // ошибка
//    {
//        printf("m=%i\tn=%i\n", m, n);
//        exit(0);
//    }
//
//    Atom3D* as = particles.GetData();
//    Compaction(as, n, area);
//}
//
////  ‘иксируем длину св€зи (принимаем текущую длину св€зи за длину недеформированной св€зи)
//void Bond3Dv2::FixLength(Atom3D* as, Area3d area, double d_av, double C_bond)
//{
//    length = float(dR(as, area).Abs());             //  длина недеформированной св€зи
//    C=float(C_bond*pow(d_av/length, GetL_pow()));   //  обновл€ем жесткость св€зи
//}
//
//
//void DetectBonds1(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, IArray<Bond3Dv2>& bonds, int& m,
//                  _List<Atom3D>& ListSlidingSealingY);
//
//void DetectBonds(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, Array<Bond3Dv2>& bonds, double d_av,
//                 int *N_layersY, double *C_bond, _List<Atom3D>& ListSlidingSealingY)
//{
//    IArray<Bond3Dv2> ibonds;    //  используем бесконечный массив
//    int m;
//    if (!area.IsPeriodic())   { print("ќбласть должна быть периодической\n"); return; }
//
//    DetectBonds1(as, n, cells, area, ibonds, m, ListSlidingSealingY);
//
//    int n_layersY0=N_layersY[0];
//    int n_layersY1=n_layersY0 + N_layersY[1];
//    int n_layersY =n_layersY1 + N_layersY[2];
//
//    bonds.Create(m);
//    int i;
//    for (i=0; i<m; i++)
//    {
//        bonds[i] = ibonds[i];
//
//        double y1 = as[bonds[i].First() ].R().y;
//        double y2 = as[bonds[i].Second()].R().y;
//        double Y = ((y1+y2)/2 - area.Y1())/area.Height()*n_layersY;
//
//        if      (Y < n_layersY0)    { bonds[i].FixLength(as, area, d_av, C_bond[0]); }
//        else if (Y < n_layersY1)    { bonds[i].FixLength(as, area, d_av, C_bond[1]); }
//        else                        { bonds[i].FixLength(as, area, d_av, C_bond[2]); }
//    }
//}
//
//void DetectBonds1(Atom3D* as, int n, ArrayS< _List<Atom3D> >& cells, Area3d area, IArray<Bond3Dv2>& bonds, int& m,
//                  _List<Atom3D>& ListSlidingSealingY)
//{
//    double		ax = area.Width(),		ay = area.Height(),		az = area.Depth();
//    double		x0 = area.X1(),			y0 = area.Y1(),			z0 = area.Z1();
//
//    Atom3D* i;
//    Atom3D* i_end = as + n;
//    long ix, iy;
//    int iz;
//
//   	double a_cut = GetAcut();
//   	double aa_cut = a_cut * a_cut;
//
//    int nx = cells.GetCount(),   ny = cells[0L].GetCount(),   nz = cells[0L][0L].GetCount();
//	double       cx = ax / nx,   cy = ay / ny,                cz = az / nz;
//
//    m = 0;
//    for (i = as; i != i_end; i++)
//    {
//        ix = long((i->R().x - x0) / cx);
//        iy = long((i->R().y - y0) / cy);
//        iz = long((i->R().z - z0) / cz);
//
//        // ищем ближайшие частицы в этой и соседних €чейках cells
//        int _jx, _jy, _jz;
//        int _ix = ix - 1, _iy = iy - 1, _iz = iz - 1;
//        int iy_ = iy + 2, iz_ = iz + 2;
//
//        // выход€щие за пределы массива cells индексы iy не рассматриваем,
//        // периодичность по оси Y не ставитс€,
//        // т.к. по оси Y ставитс€ скольз€ща€ заделка
//        _iy = _iy < 0 ? 0 : _iy;
//        iy_ = iy_ >= ny ? ny : iy_;
//        if (iy==0 || iy==ny-1)  { ListSlidingSealingY.Insert(i); }
//
//		for(_jx = _ix;	_jx <= ix;	_jx++) // не рассматриваем соседние €чейки cells справа (исключаем дублирование св€зей - частиц i,j и j,i)
//		for(_jy = _iy;	_jy < iy_;	_jy++)
//		for(_jz = _iz;	_jz < iz_;	_jz++)
//		{
//            char shift = 0;
//            long jx = _jx, jy = _jy;
//            int jz = _jz;
//
//            // ставим периодические граничные услови€ на гран€х, перпендикул€рных ос€м X, Z
//            if (jx < 0)    { jx = nx - 1;  shift = 1; }
//            // else if (jx == nx)  { jx = 0;  shift = 2; }   (вообще jx <= ix < nx)
//            if (jz < 0)    { jz = nz - 1;  shift |= 16; }
//            else if (jz == nz)  { jz = 0;  shift |= 32; }
//
//			_List<Atom3D>* lj = &cells[jx][jy][jz];
//			BOOL more = lj->Restore();
//			while(more)
//			{
//				Atom3D* j = lj->Iterate(more);
//
//                if (i == j) { _jy = iy_; _jz = iz_; break; }    // выходим из трех циклов (по _jx, _jy, _jz) -
//                                                                // исключаем дублирование св€зей - частиц i,j и j,i
//                Bond3Dv2 bond(int(i-as), int(j-as), shift);
//                double dd = bond.dR(as, area).Sqr();
//
//                if (dd < aa_cut)
//                {
//                    bonds[m++] = bond;
//                    i->N_Bonds()++;
//                    j->N_Bonds()++;
//                }
//            }
//        }
//    }
//}
//
//void AngleBond3D::FixAngle(Bond3Dv2* bonds, double d_av, double C_bond)
//{
//    Bond3Dv2& b1 = bonds[First()];
//    Bond3Dv2& b2 = bonds[Second()];
//    Vector3D& dR1 = b1.dR();
//    Vector3D& dR2 = b2.dR();
//    double cos_fi =
//        (b2.First()==b1.First() || b2.Second()==b1.Second()) ?
//        (dR1*dR2) / dR1.Abs() / dR2.Abs():
//       -(dR1*dR2) / dR1.Abs() / dR2.Abs();
//
//    angle = float(acos(cos_fi));
//    C=float(C_bond*GetC_fi()*min2(dR1.Sqr(), dR2.Sqr()) / d_av / d_av);  // alpha=2
//}
//
//void DetectAngleBonds1(int n, Bond3Dv2* bonds, ArrayS<int>& group,
//                      Array<AngleBond3D>& a_bonds);
//
//void DetectAngleBonds(Atom3D* as, int n, Bond3Dv2* bonds, ArrayS<int>& group, Area3d area,
//                      Array<AngleBond3D>& a_bonds, double d_av, int *N_layersY, double *C_bond)
//{
//    DetectAngleBonds1(n, bonds, group, a_bonds);
//
//    int n_layersY0=N_layersY[0];
//    int n_layersY1=n_layersY0 + N_layersY[1];
//    int n_layersY =n_layersY1 + N_layersY[2];
//
//    int j, na_bonds = a_bonds.GetCount();
//    for (j = 0; j < na_bonds; j++)
//    {
//        int i1 = a_bonds[j].First();
//        int i2 = a_bonds[j].Second();
//        int i3 = bonds[i1].First();
//        int i4 = bonds[i1].Second();
//        int i5 = bonds[i2].First();
//        int i6 = bonds[i2].Second();
//        int i7 = (i3==i5 || i3==i6) ? i3 : i4;
//        double Y = (as[i7].R().y - area.Y1())/area.Height()*n_layersY;
//
//        if      (Y < n_layersY0)    { a_bonds[j].FixAngle(bonds, d_av, C_bond[0]); }
//        else if (Y < n_layersY1)    { a_bonds[j].FixAngle(bonds, d_av, C_bond[1]); }
//        else                        { a_bonds[j].FixAngle(bonds, d_av, C_bond[2]); }
//    }
//}
//
//// три сло€ по оси Y с разными упругими модул€ми
//void CreateConfig3D(Array<Atom3D>& particles, int n_layersX, int *N_layersY, int n_layersZ,
//                    double* C_bond, Area3D& area, double d_av,
//                    Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, _List<Atom3D>& ListSlidingSealingY)
//{
//    int n_layersY=N_layersY[0]+N_layersY[1]+N_layersY[2];
//    //  массив частиц со случайными отклонени€ми от плотноупакованной решетки
//    NewSystem3D_01(particles, n_layersX, n_layersY, n_layersZ, area, d_av);
//
//    Atom3D* as = particles.GetData();
//    int n = particles.GetCount();
//    ArrayS< _List<Atom3D> > cells;
//
//    //  сортировка частиц по €чейкам
//    CreateMap(as, n, cells, area);
//
//    //  нахождение массива св€зей
//    DetectBonds(as, n, cells, area, bonds, d_av, N_layersY, C_bond, ListSlidingSealingY);
//
//    if (GetC_fi()!=0.0)
//    {
//        //  находим группы св€зей при каждой частице
//        ArrayS<int> group;
//        DetectGroup(as, n, bonds, group);
//
//        //  нахождение массива углов между св€з€ми
//        DetectAngleBonds(as, n, bonds.GetData(), group, area, a_bonds, d_av, N_layersY, C_bond);
//    }
//}

