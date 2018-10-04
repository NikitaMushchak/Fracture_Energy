// ������ ������
#include "atom3Dv2.h"
#include "util.h"
#include "a3r.h"
#include "a3r_b.h"
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>

void Modeling3D_10_0();
void Modeling3D_10_1();
void Modeling3D_10_2();
void Modeling3D_11_1();
void Modeling3D_14();
void Modeling3D_14eps();
void Modeling3D_15eps();
void Modeling3D_17();
void Modeling3D_16eps();
void Modeling3D_seq();
void Modeling3D_18();

void Modeling3Dv2()
{
    //Modeling3D_10_0();    // ���� - �������� ��������� ������������ � ����� �������� a3r
    //Modeling3D_10_1();    // ���������� ������� ������� (�����), ��� ����������
    //Modeling3D_10_2();    // ���������� ������� �������, ���������� (����������)
    //Modeling3D_11_1();    // ���������� ��������� �����������, ���������� (����������), �����
    Modeling3D_14();      // 
    //Modeling3D_14eps();   // ����������� �������� ��� ������ ����������
    //Modeling3D_15eps();   // �� ��, ��� �������
    //Modeling3D_17();      // ������� - ���������� (���������)
    //Modeling3D_16eps();   // ����������, ���� ������
    //Modeling3D_seq();     // ���������������� ���� ���������� ������
    //Modeling3D_18();      // ��������� ����� ������
}

static double d_av = 1; // ������� ���������� ����� ���������

// ���� - �������� ��������� ������������ � ����� �������� a3r
void Modeling3D_10_0()
{
    int random_init=(int)time(0);

    SetAngleSprings(0);
    srand(random_init);

    Array<Atom3D> particles;
    int n_layers = 20;
    Area3D area;

    SetTypeMaterial(1);

    area.Periodic() = 1;            // ������������� ��������� �������
    Array<Bond3Dv2> bonds;          // ������ ������ ����� ���������
    Array<AngleBond3D> a_bonds;     // ������ ����� ����� �������

    CreateConfig3D(particles, n_layers, area, d_av, bonds, a_bonds);
    //CreateConfig3D(particles, n_layers, area, d_av, bonds, a_bonds, 0.3);

    for (int i=0; i<particles.GetCount(); i++)  { particles[i].Color()=(char)(i%3); }
    char str[100] = "output/a0.a3r";
    Save_A3R_B(str, particles.GetData(), particles.GetCount(), d_av/2, 0, 1);

    printf("Sizeof Atom3D=%i\n", sizeof(Atom3D));
    printf("Sizeof Bond3Dv2=%i\n", sizeof(Bond3Dv2));
    printf("Sizeof AngleBond3D=%i\n", sizeof(AngleBond3D));
    return;
}

// ���������� ���������� ������
double DefBond(Atom3D* as, Array<Bond3Dv2>& bonds, const Area3D& area, double scal, Array<Atom3D>& points0)
{
    int N_B = bonds.GetCount();

    double eps = 0;
    int i;
    for (i = 0; i < N_B; i++)
    {
        double L0 = bonds[i].Length();
        double L1 = bonds[i].dR(points0.GetData(), area).Abs();
        double L2 = bonds[i].dR(as, area).Abs();

        eps += sqr((L2-L1)/L0);
    }

    eps /= N_B;
    eps = sqrt(eps);
    eps /= (scal-1);

    print("����. ����������: %f\n", eps);
    return eps;
}

static int outp=0;
static int n_layers = 20;//40 // ������ ������� (� ������� ����� ������)

// ���������� ������� �������

double Modeling3D_10(double* pC11=0, double* pC12=0, double* p2C44=0, double* pNu=0)
{
    Array<Atom3D> particles;
    Area3D area;

  //SetTypeMaterial(2, 1-rnd_dspl);
  //SetTypeMaterial(3, 1-rnd_dspl); ���������, ���������� �������� Modeling3D_11_1()

    //srand((int)time(0));          // ����� ���������������� ����� ������� Modeling3D_10!

    area.Periodic() = 1;            // ������������� ��������� �������
    Array<Bond3Dv2> bonds;          // ������ ������ ����� ���������
    Array<AngleBond3D> a_bonds;     // ������ ����� ����� �������

    CreateConfig3D(particles, n_layers, area, d_av, bonds, a_bonds);
    //CreateConfig3D(particles, n_layers, area, d_av, bonds, a_bonds, 0.5); // ������ (0.5 ���. ~ 30 ���.)

    Bond3Dv2* bonds_d = bonds.GetData();
    int N_B = bonds.GetCount();
    Atom3D* as = particles.GetData();
    int n = particles.GetCount();
    AngleBond3D* a_bonds_d = a_bonds.GetData();
    int N_AB = a_bonds.GetCount();

    print("���������� ������: %i\n", n);
    print("������/������: %f\n", (double)N_B/n);
    print("�����/������: %f\n", (double)N_AB/n);

    //double scal = 1.01;
    double scal = 1.0001;

    // ��������� ��������� ������
    Array<Vector3D> points_b(n);
    int i;
    for (i=0; i<n; i++) { 
		points_b[i] = as[i].R(); }
    Area3D area0 = area;

    // ����������� ������ ����� ��� X
    area.ExtendX(as, n, scal);

    // ��������� ��������� ������ ����� ����������
    Array<Atom3D> points0(n);
    for (i=0; i<n; i++) { points0[i] = as[i]; }

    //int outp = 0;
    FILE* f;
    //if (outp) { f = fopen("���.��.txt", "w"); }

    // ���������� ����������

    SetB_dt(1.0);
    int s=0;
    do
    {
		std::cout << "Relaxing\n";
        //Step3D(as, n, area, bonds_d, N_B, a_bonds_d, N_AB);
        if (outp)
        {
            //fprintf(f, "%i\t%e\t", s, GetV2(as, n));
            //fprintf(f, "%e\t", GetLastAvDisplacement(as, n));
            ////fprintf(f, "%e\n", GetFullAvDisplacement(as, n, points0.GetData())/3/s);
            //fprintf(f, "%e\n", GetFullAvDisplacement(as, n, points0.GetData()));
        }
        if (!(s%100))
        {
            print("���: %i\n", s);
            DefBond(as, bonds, area, scal, points0);
        }
        s++;

    } while(GetLastAvDisplacement(as, n)>GetFullAvDisplacement(as, n, points0.GetData())/3/s && s<100);

//    if (outp) { fclose(f); }

    // ���������� ���������� ������

    double eps = 0;
    for (i = 0; i < N_B; i++)
    {
        double L0 = bonds[i].Length();
        double L1 = bonds[i].dR(points0.GetData(), area).Abs();
        double L2 = bonds[i].dR(as, area).Abs();

        eps += sqr((L2-L1)/L0);
    }

    eps /= N_B;
    eps = sqrt(eps);
    eps /= (scal-1);

    print("����. ����������: %f\n", eps);

    Tensor tau;

    // ���������� ����������
    tau = GetStressTensor(as, n, area, bonds.GetData(), N_B, tau.x.x, tau.y.y, tau.z.z, tau.x.y, tau.x.z, tau.y.z);

	std::cout << "Stresses" << std::endl;
	printf ("\n\n");
    printf ("%f\t%f\t%f\n", tau.x.x, tau.x.y, tau.x.z);
    printf ("%f\t%f\t%f\n", tau.y.x, tau.y.y, tau.y.z);
    printf ("%f\t%f\t%f\n", tau.z.x, tau.z.y, tau.z.z);/**/

    double C11 = tau.x.x/(scal-1);
    double C12 = (tau.y.y+tau.z.z)/(scal-1)/2;
    double C11_C12 = C11-C12;
    print ("\n\n");
    print ("C11= %f\n", C11);
	std::cout << "C11=" << C11 << std::endl;
	std::cout << "C12=" << C12 << std::endl;
    print ("C12= %f\n", C12);
    print ("C11-C12= %f\n", C11_C12);
    print ("����. ��������: %f\n", C12/(C11+C12));
	std::cout << "Poisson=" << C12 / (C11 + C12) << std::endl;
    if (pC11) { *pC11=C11; }
    if (pC12) { *pC12=C12; }
    if (pNu)  { *pNu=C12/(C11+C12); }

    //system("pause");

    // ��������� ��������� ������ ����� �������
    for (i=0; i<n; i++) { as[i].R() = points_b[i]; as[i].V() = VECT3D_NULL; }

    // ��������������� �������������� ������
    area = area0;
    //double aw = area.Width();

    //double ShearYX = 0.01;
    double ShearYX = 0.0001;   // ShearYX == 2*Epsilon_YX

    // �������� ������ � ����������� �� Y ����� ��� X
    area.AddShearYX(as, n, ShearYX);

    // ��������� ��������� ������ ����� ������
    for (i=0; i<n; i++) { points0[i] = as[i]; }

    //if (outp) { f = fopen("���.��.2.txt", "w"); }

    // ���������� ����������

    s = 0;
    do
    {
        //Step3D(as, n, area, bonds_d, N_B, a_bonds_d, N_AB);
        //s++;
        if (outp)
        {
            //fprintf(f, "%i\t%e\t", s, GetV2(as, n));
            //fprintf(f, "%e\t", GetLastAvDisplacement(as, n));
            ////fprintf(f, "%e\n", GetFullAvDisplacement(as, n, points0.GetData())/3/s);
            //fprintf(f, "%e\n", GetFullAvDisplacement(as, n, points0.GetData()));
        }
        if (!(s%100)) { print("���: %i\n", s); }
        s++;

    } while(GetLastAvDisplacement(as, n)>GetFullAvDisplacement(as, n, points0.GetData())/3/s && s<100);

    //if (outp) { fclose(f); }

    // ���������� ����������
    //tau = GetStressTensor(as, n, area, bonds.GetData(), N_B);

    /*printf ("\n\n");
    printf ("%f\t%f\t%f\n", tau.x.x, tau.x.y, tau.x.z);
    printf ("%f\t%f\t%f\n", tau.y.x, tau.y.y, tau.y.z);
    printf ("%f\t%f\t%f\n", tau.z.x, tau.z.y, tau.z.z);*/

    double _2C44 = (tau.x.y+tau.y.x)/ShearYX;
    print("\n\n");
    std::cout << "2*C44= %f\n" << _2C44;
    print("�������� �����������: %f\n", _2C44/C11_C12);

    if (p2C44) { *p2C44=_2C44; }

    return _2C44/C11_C12;
}
/**/
//
//// ���������� ������� ������� (�����), ��� ����������
void Modeling3D_10_1()
{
    FILE* f = fopen("������1.txt", "a");

    //double L_pow = 2.82; // ��� d_min==1
    //double L_pow = 2.8425; // ��� d_min==0.9
    //double L_pow = 1.85; //- ��� C_fi == -0.07
    //double L_pow = 1; // ��� type_material==4
    //double A_cut = 1.9;
    double C11;
    double C12;
    double _2C44;
    double Nu;
    double an;
    int random_init=(int)time(0);

	n_layers = 30;

    int i;
  //for (i=0; i<10; i++)
    for (i=-1; i<=-1; i++)
    {
        //SetTypeMaterial(3, 0.9999999, &L_pow);
        //SetTypeMaterial(3, 0.9, &L_pow);
        //SetTypeMaterial(3, 0.82, &L_pow);
        //SetTypeMaterial(4, 0.394, &L_pow);
        //SetTypeMaterial(2, 0.33, &L_pow, &A_cut);
        SetTypeMaterial(0.33);

        //double C_fi = -0.11;
        double C_fi = 0;
        SetAngleSprings(C_fi);
        srand(random_init);
        an = Modeling3D_10(&C11, &C12, &_2C44, &Nu);
        printf    ("%f\t%f\t%f\t%f\t%f\t%f\n", C_fi, C11, C12, _2C44, Nu, an);
        fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\n", C_fi, C11, C12, _2C44, Nu, an);
    }
    fclose(f);
}

// ���������� ������� �������, ���������� (����������)
void Modeling3D_10_2()
{
    double C11[10], C11_av=0, C11_2=0;
    double C12[10], C12_av=0, C12_2=0;
    double _2C44[10], _2C44_av=0, _2C44_2=0;
    double Nu[10], Nu_av=0, Nu_2=0;
    double an[10], an_av=0, an2=0;

    //double L_pow = 1;
    //double A_cut = 1.9;

  SetTypeMaterial(0.33);
    //SetTypeMaterial(0.999);
    SetAngleSprings(0);
  //SetAngleSprings(4);

    int j, n = 10;
    for (j=0; j<n; j++)
    {
        //if (j==0) { outp=1; } else { outp=0; }
        srand((int)time(0));
        an[j] = Modeling3D_10(C11+j, C12+j, _2C44+j, Nu+j);  // �������� 2
        //an[j] = Modeling3D_10(0.18, C11+j, C12+j, _2C44+j, Nu+j);  // �������� 3
        an_av += an[j];
        C11_av += C11[j];
        C12_av += C12[j];
        _2C44_av += _2C44[j];
        Nu_av += Nu[j];
    }
    an_av/=n;
    C11_av/=n;
    C12_av/=n;
    _2C44_av/=n;
    Nu_av/=n;
    for (j=0; j<n; j++)
    {
        an2+=sqr(an[j]-an_av);
        C11_2+=sqr(C11[j]-C11_av);
        C12_2+=sqr(C12[j]-C12_av);
        _2C44_2+=sqr(_2C44[j]-_2C44_av);
        Nu_2 +=sqr(Nu[j]-Nu_av);
    }
    an2/=n;   an2=sqrt(an2)*M_SQRT2;
    C11_2/=n;   C11_2=sqrt(C11_2)*M_SQRT2;
    C12_2/=n;   C12_2=sqrt(C12_2)*M_SQRT2;
    _2C44_2/=n;   _2C44_2=sqrt(_2C44_2)*M_SQRT2;
    Nu_2/=n;    Nu_2= sqrt(Nu_2)*M_SQRT2;
    //double an_av=Modeling3D_1(0.7);
    printf("an = %f  +/- %f\n", an_av, an2);
    printf("C11= %f  +/- %f\n", C11_av, C11_2);
    printf("C12= %f  +/- %f\n", C12_av, C12_2);
    printf("2C44= %f  +/- %f\n", _2C44_av, _2C44_2);
    printf("Nu= %f  +/- %f\n", Nu_av, Nu_2);
}

//// ���������� ��������� �����������, ���������� (����������), �����
//void Modeling3D_11_1()
//{
//    FILE* f;
//    //int i;
//    double d_min_d = 0.999;//0.35;
//    //for (i=0; i<19; i++)
//    for (; d_min_d<1.01;)
//    {
//        double an[40], an_av=0, an2=0;
//      //double d_min_d = 0.33;//0.384;//0.394;//0.21;
//      //double L_pow = 2.82 + 0.01*i;
//      //double L_pow = 2.83 + 0.01*i/3;
//      //double L_pow = 1;
//      //double A_cut = 1.9;
//
//        //SetTypeMaterial(d_min_d, L_pow, A_cut);
//        SetTypeMaterial(d_min_d);
//
//        //int j, n = i==0 ? 1 : 10;
//        int j, n = 10;
//        //int j, n = 40;
//        for (j=0; j<n; j++)
//        {
//            srand((int)time(0));
//            an[j] = Modeling3D_10();
//            an_av += an[j];
//        }
//        an_av/=n;
//        for (j=0; j<n; j++) { an2+=sqr(an[j]-an_av); }
//        an2/=n;
//        an2=sqrt(an2);
//
//        f = fopen("�����.����.2.txt", "a");
//      //printf("\t\t\t%f\t%f +/- %f\n", L_pow, an_av, an2);
//      //fprintf(f, "%f\t%f\t%f\n", L_pow, an_av, an2);
//        printf("\t\t\t%f\t%f +/- %f\n", d_min_d, an_av, an2);
//        fprintf(f, "%f\t%f\t%f\t%f\n", d_min_d, 1-d_min_d, an_av, an2);
//        fclose(f);
//
//        d_min_d += 0.05;
//    }
//}

/////////////////////////////////////////////////////////////////////////////////////
// ��������� �������

// ���������� ����� ���������� ������ ����� ��� X
/*
void Equlibr1(Array<Atom3D>& particles, Area3D& area, Array<Bond3Dv2>& bonds,
              Array<AngleBond3D>& a_bonds, int s, Array<Vector3D>& points, double V2_max)
{
    Atom3D* as = particles.GetData();
    int n = particles.GetCount();
    Bond3Dv2* bonds_d = bonds.GetData();
    int N_B = bonds.GetCount();
    AngleBond3D* a_bonds_d = a_bonds.GetData();
    int N_AB = a_bonds.GetCount();

    double V2=GetV2(as, n);        // ����� ��������� �������� ��������
    SetB_dt(1.0);                  // ������

    for (; V2*1e4>=V2_max; s++)
    {
        if (!(s%500))
        {
            SaveConfig("cnf1.tmp", particles, area, bonds, a_bonds, s, &V2_max);
            system("del cnf1.dat");
            system("ren cnf1.tmp cnf1.dat");
        }

        if (!(s%100))
        {
            print("���: %i\n", s);
            //print("������: %i\n", bonds.GetCount()+bonds_liq.GetCount());
            print("���.�������: %e\n", V2);
            print("���.������� (����): %e\n", V2_max);
            print("����� �������: %f\n", Volume(bonds));
            print("�����: %.1f\n", clock()/1000.);
        }

        if (!(s%10))
        {
            FILE* f=fopen("output/���������_����������1.txt", "a");
            fprintf(f, "%i\t%e\t%e\n", s, V2, GetFullAvDisplacement(as, n, points.GetData()));
            fclose(f);
        }

        Step3D(as, n, area, bonds_d, N_B, a_bonds_d, N_AB);
        V2=GetV2(as, n);
        V2_max=max2(V2, V2_max);
    }
}
/**/
////////////////////////////////////////////////////
// ��������� �������

// ���������� ����� ���������� �������
/*
void Equlibr2(Array<Atom3D>& particles, Area3D& area, Array<Bond3Dv2>& bonds,
              Array<AngleBond3D>& a_bonds, int s, Array<Vector3D>& points, double V2_max)
{
    Atom3D* as = particles.GetData();
    int n = particles.GetCount();
    Bond3Dv2* bonds_d = bonds.GetData();
    int N_B = bonds.GetCount();
    AngleBond3D* a_bonds_d = a_bonds.GetData();
    int N_AB = a_bonds.GetCount();

    double V2=GetV2(as, n);        // ����� ��������� �������� ��������
    SetB_dt(0.1);                  // ������

    for (; V2*1e3>=V2_max; s++)
    //for (; V2*1e2>=V2_max; s++)
    {
        if (s && !(s%500))
        {
            SaveConfig("cnf1.tmp", particles, area, bonds, a_bonds, s, &V2_max);
            system("del cnf1.dat");
            system("ren cnf1.tmp cnf1.dat");
        }

        if (!(s%100))
        {
            print("���: %i\n", s);
            //print("������: %i\n", bonds.GetCount()+bonds_liq.GetCount());
            print("���.�������: %e\n", V2);
            print("���.������� (����): %e\n", V2_max);
            print("����� �������: %f\n", Volume(bonds));
            print("�����: %.1f\n", clock()/1000.);

            //GetStress(as, n, area, bonds, 1, 0.25).x;
        }

        if (!(s%10))
        {
            FILE* f=fopen("output/���������_����������2.txt", "a");
            fprintf(f, "%i\t%e\t%e\n", s, V2, GetFullAvDisplacementGap(as, n, points.GetData()));
            fclose(f);
        }

        Step3D(as, n, area, bonds_d, N_B, a_bonds_d, N_AB);
        V2=GetV2(as, n);
        V2_max=max2(V2, V2_max);
    }
}
/**/
//double Modeling3D_16(double size_cr, int n_layers, int continued)
//{
//    Array<Atom3D> particles;
//    Area3D area;
//
//  //SetTypeMaterial(1, 0.22);
//    SetTypeMaterial(0.30);
//  //SetTypeMaterial(3, 0.82);
//
//    //time_t tm;
//    //printf("time()=%i\n", (int)time(&tm));
//    srand((int)time(0));
//
//    area.Periodic() = 1;            // ������������� ��������� �������
//    Array<Bond3Dv2> bonds;          // ������ ������ ����� ���������
//    Array<AngleBond3D> a_bonds;     // ������ ����� ����� �������
//
//    if (!continued)  { CreateConfig3D(particles, n_layers, area, d_av, bonds, a_bonds); }
//
//    // ����������� ������ ����� ��� X
//    double scal=1.0001;//1.00001
//    if (!continued)  { area.ExtendX(particles.GetData(), particles.GetCount(), scal); }
//
//    int s=0;
//    double V2_max=0;               // ����� ��������� �������� �������� (��������)
//
//    if (continued)   { LoadConfig("cnf1.dat", particles, area, bonds, a_bonds, s, &V2_max); }
//
//    Atom3D* as = particles.GetData();
//    int n = particles.GetCount();
//
//    int N_equ=1;
//    // ��������� ��������� ������ ����� ����������
//    Array<Vector3D> points0;
//    if (!continued)
//    {
//        points0.Create(n);
//        for (int i=0; i<n; i++) { points0[i] = as[i].R(); }
//        SaveConfig2("cnf2.dat", points0, 1);
//    }
//    else { LoadConfig2("cnf2.dat", points0, N_equ); } // N_equ ����� ���� == 1 ��� 2
//
//    FILE* f;
//    if (!continued)  { f=fopen("output/���������_����������1.txt", "w"); fclose(f); }
//
//    // ������� ���������� ����� ���������� ������ ����� ��� X
//    if (N_equ==1)  { Equlibr1(particles, area, bonds, a_bonds, s, points0, V2_max); }
//
//    Bond3Dv2* bonds_d = bonds.GetData();
//    int N_B = bonds.GetCount();
//
//    double Sigm0=0;
//
//    if (N_equ==1) // ������ ���� ��� �� ���������
//    {
//        Sigm0=GetStress(as, n, area, bonds_d, N_B).x.x;
//
//        Vector3D center(0.5,0.5,0.5);
//        Vector3D normal(1,0,0);
//
//        // �������� �������
//        CreateFracture(as, area, bonds, center, size_cr, normal);
//
//        Bond3Dv2* pB = bonds_d;
//        Bond3Dv2* pE = pB+N_B;
//        const char thr=Bond3Dv2::GapM();
//        for (pB=bonds_d; pB<pE; pB++)      // ��� �� ����������� ����� - �������������
//        {
//            if (pB->Gap() < thr) { pB->Gap()=0; }
//        }
//
//        print("\n�������� �������\n\n");
//
//        // ��������� ���������� ������ �� ���� ������� ����� ����������� ����������
//        Save_A3R_Gap("output/���������.a3r", as, n, bonds, area, 0.5);
//
//        V2_max=0; s=0;
//
//        for (int i=0; i<n; i++) { points0[i] = as[i].R(); }
//        SaveConfig("cnf1.tmp", particles, area, bonds, a_bonds, s, &V2_max);
//        SaveConfig2("cnf2.tmp", points0, 2);
//        system("del cnf1.dat");
//        system("del cnf2.dat");
//        system("ren cnf1.tmp cnf1.dat");
//        system("ren cnf2.tmp cnf2.dat");
//
//        f=fopen("output/���������_����������2.txt", "w");
//        fclose(f);
//    }
//
//    // ������� ���������� ����� ���������� ������� (������ ����)
//    Equlibr2(particles, area, bonds, a_bonds, s, points0, V2_max);
//
//    Save_A3R_Gap("output/���.a3r", as, n, bonds, area, 0.5);
//
//    //double Sigm1=GetStress(as, n, area, bonds, 1, 0.25).x;
//    double Sigm1=GetStress(as, n, area, bonds_d, N_B).x.x;
//    print("����������X ��� �������: %e\n", Sigm0);
//    print("����������X: %e\n", Sigm1);
//
//    return Sigm1;
//}

//#include "����.h"
//
//// ������� - ���������� (���������)
//void Modeling3D_17()
//{
//    double size_cr = 0.2;//0.35;    // ������������� ������ �������
//    int    n_layers= 80;     //50
//    int continued = 1;				// ����������� ���������� � �������������� ������ �������� ��������� cnf1.dat � cnf2.dat
//
//    double Sigm1 = Modeling3D_16(size_cr, n_layers, continued);
//    int n, i;
//    double r;
//    Vector3D* as0=Load_A3R("output/���������.a3r", n, r);
//    Vector3D* as =Load_A3R("output/���.a3r", n, r);
//    Array<double> bn_e(n), R(n);
//    double a=(n_layers*M_SQRT2)*size_cr/2;
//  //double bn_coef=4/M_PI*a/0.982*1.113e-4; // 0.982=E/(1-nu^2)=(C11-C12)/(1-nu) ��� ��������� 2
//  //                                        // 1.113e-4=p
//    double bn_coef=4/M_PI*a/0.982*Sigm1;    // 0.982=E/(1-nu^2)=(C11-C12)/(1-nu) ��� ��������� 2
//
//    for (i=0; i<n; i++)
//    {
//        R[i]=sqrt(sqr(as[i].z)+sqr(as[i].y))/a;
//        bn_e[i]=fabs(as[i].x-as0[i].x)/bn_coef;
//    }
//
//    Array<int> ind(n);
//    SortInit(ind.GetData(), n);
//    SortSplit(ind.GetData(), n, R.GetData());
//
//    FILE* f=fopen("output/���������.txt", "w");
//    for (i=0; i<n; i++)
//    {
//        double R_=R[ind[i]];
//        fprintf(f, "%f\t%f\t%f\n", R_, bn_e[ind[i]], R_>=1? 0 : sqrt(1-R_*R_));
//        //fprintf(f, "%f\t%f\n", as[ind[i]].x/bn_coef, as0[ind[i]].x/bn_coef); - ������ Save_A3R_B(file_name, prt_bond, r, 0, 1);
//    }
//    fclose(f);
//
//    delete [] as0;
//    delete [] as;
//}

// ��������� ����� ������
//void Modeling3D_18()
//{
//    Array<Atom3D> particles;
//    int n_layers = 40;//10
//    Area3D area;
//    //double L_pow = 1;
//    //double A_cut = 1.9;
//
//    SetTypeMaterial(0.33); // L_pow = 1;
//    SetAngleSprings(0);
//
//    srand((int)time(0));
//
//    area.Periodic() = 1;            // ������������� ��������� �������
//    Array<Bond3Dv2> bonds;          // ������ ������ ����� ���������
//    Array<AngleBond3D> a_bonds;     // ������ ����� ����� �������
//
//    CreateConfig3D(particles, n_layers, area, d_av, bonds, a_bonds);//, 0.4);//, 0.78);
//
//    Atom3D* as = particles.GetData();
//    int n = particles.GetCount();
//    int N_B = bonds.GetCount();
//    //int N_AB = a_bonds.GetCount();
//
//    double nx[20], nx_=0, nx2=0;
//    int i;
//    for (i=0; i<20; i++)
//    {
//        nx[i] = GetNBonds(as, n, area, bonds.GetData(), N_B, 1, 0.5+(i-10)/100.0);
//        nx_ += nx[i];
//    }
//    nx_ /= 20;
//
//    for (i=0; i<20; i++)
//    {
//        nx2 += sqr(nx[i]-nx_);
//    }
//    nx2 /= 20;
//    nx2 = sqrt(nx2);
//
//    print("������*cos(alpha) �� ������� �������: %f\n", nx_);
//    print("������� ����������: %f\n", nx2);
//    print("������: %i\n", N_B);
//}

void UpdateParams(int& s_max, double& press_k)
{
    FILE* f = fopen("���������1.txt", "r");
    if(f==0) { print("����������� ���� \"���������1.txt\"\n"); s_max=0; return; }
    int err=0;
    s_max=InputInt(f, err);
    press_k=InputDouble(f, err);
	fclose(f);

    if(err) { print("������ UpdateParams\n"); }
}

void UpdateParams_E(int& s_max, double& press_k, double& mod_E, double& KC, int& n_layers, 
	double& part_center, double& diametr)
{
	/*FILE* f = fopen("���������1.txt", "r");
	if (f == 0) { print("����������� ���� \"���������1.txt\"\n"); s_max = 0; return; }
	int err = 0;
	s_max = InputInt(f, err);
	press_k = InputDouble(f, err);
	mod_E = InputDouble(f, err);
	KC = InputDouble(f, err);
	n_layers = InputInt(f, err);
	part_center = InputDouble(f, err);
	diametr = InputDouble(f, err);
	fclose(f);*/

	//if (err) { print("������ UpdateParams\n"); }
	FILE *file;
	file = fopen("���������1.txt", "r");
	//string pr;
	char str[80];

	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", s_max);
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", press_k);
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", mod_E);
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", KC);
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", n_layers);
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", part_center);
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", diametr);
	fclose(file);
	std::cout << s_max << " " << press_k << " " << mod_E << " " << KC << " " << n_layers << " " << part_center << " " << diametr << "\n";
}

void CreateFracture_1(Atom3D* as, Area3d area, Array<Bond3Dv2>& bonds, Vector3D center, double diameter, Vector3D normal)
{
    
	/*if (orientation == 1){
		normal.x = 1;
	}else if (orientation == 2){
		normal.z = 1;
		normal.x = 1;
	}else if(orientation == 3){
		normal.z = 1;
	}*/

    // �������� �������
    CreateFracture(as, area, bonds, center, diameter, normal);
}

void CreateFracture_2(Atom3D* as, Area3d area, Array<Bond3Dv2>& bonds)
{
    Vector3D center(0.5,0.5,0.4);
    Vector3D normal(0,0,1);

    // �������� ������
    CreateFracture(as, area, bonds, center, 0.25, normal);
    center.Set(0.5,0.5,0.6);
    CreateFracture(as, area, bonds, center, 0.25, normal);
}

void CreateFracture_mass(Atom3D* as, Area3d area, Array<Bond3Dv2>& bonds, double density, double D, int parallel, int intersection,int orientation)
{
	Vector3D center(0.5, 0.5, 0.5);
	Vector3D normal(0, 0, 0);
	if (orientation == 1){
		normal.x = 1;
	}else if (orientation == 2){
		normal.z = 1;
		normal.x = 1;
	}else if (orientation == 3){
		normal.z = 1;
	}
	Vector3D r(0, 0, 0); // ������, ����������� ������ ������

	double cos1; // ������� ���� ����� �������� ������ ������� � r
	double cos2; // ������� ���� ����� �������� dnjhjq ������� � r

	double S = D * D * D / 8;
	int N = ceil(density * 1 / S); // ���������� ������
	int k = 0; // ������� ���������� ������ ��� �����
	int check = 1; // �������� �� �������������

	print("���������� ������: %i\n", N);

	std::vector<Vector3D> coord(N);
	std::vector<Vector3D> fracnorm(N);

	// �������� ������
	for (int i = 1; i <= N; i++)
	{
		//print("K: %i\n", k);
		check = 1;

		//������������� �������
		if (parallel == 0)
		{
			normal.Set((double)(rand()) / RAND_MAX*(1 - (-1)) + (-1), (double)(rand()) / RAND_MAX*(1 - (-1)) + (-1), (double)(rand()) / RAND_MAX*(1 - (-1)) + (-1));
		}

		//������������� �����
		center.Set((double)(rand()) / RAND_MAX*(1 - D / 2 - (0 + D / 2)) + (0 + D / 2), (double)(rand()) / RAND_MAX*(1 - D / 2 - (0 + D / 2)) + (0 + D / 2), (double)(rand()) / RAND_MAX*(1 - D / 2 - (0 + D / 2)) + (0 + D / 2));

		//������� �� ���������
		if (intersection == 0)
		{
			for (int j = 0; j <= k; j++)
			{
				r = coord[j] - center;

				if (r.Abs() < D)
				{
					cos1 = r * normal / (r.Abs() * normal.Abs());
					cos2 = r * fracnorm[j] / (r.Abs() * fracnorm[j].Abs());

					if (abs(D * sqrt(1 - cos1 * cos1) / 2 + D * sqrt(1 - cos2 * cos2) / 2) + 0 > r.Abs())
					{
						check = 0;
					}
				}

			}

		}

		if (check == 1)
		{
			CreateFracture(as, area, bonds, center, D, normal);
			coord[k] = center;
			fracnorm[k] = normal;
			k += 1;
		}
		else
		{
			N = N + 1;
		}
	}
}

// 
void Modeling3D_14()
{
	/*=========================*/
	int s = 0, s0 = 0, outp = 1, s_max = 30000;
	double press = 0;
	double press_k = 0;
	double press_k_integ = 0;
	double mod_E_in = 1;
	double KC_init = 1;
	double mod_E = 1;
	double KC = 1;
	double part_third_centr = 0.8;
	Array<Atom3D> particles;
	int n_layers = 0;//10
	Area3D area;
	double diameter = 0.2;
	
	double sigma_x = 0;
	double sigma_y = 0;
	double sigma_z = 0;
	double sigma_xy = 0;
	double sigma_xz = 0;
	double sigma_yz = 0;
	double RealP = 0;
	int Number_crack = 0;
	double orientation = 0;
	int parallel = 0;
	int tt = 0;
	size_t v = 5; // compression step v * sexit

	FILE* f1 = fopen("output/Crack_Propagation.txt", "w");
	fclose(f1);
	FILE* f2 = fopen("output/Stresses.txt", "w");
	fclose(f2);
	FILE* f3 = fopen("output/Maximum.txt", "w");
	fclose(f3);
	FILE* f4 = fopen("output/Real_Press.txt", "w");
	fclose(f4);
	FILE* f5 = fopen("output/Real_Press_L.txt", "w");
	fclose(f5);

	//UpdateParams_E(s_max, press_k, mod_E, KC,n_layers, part_third_centr,diameter); // press_k=1e-3;
	FILE *file;
	file = fopen("data.dat", "r");
	
	// // //string pr;
	char str[80];

	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", &(s_max));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", &(press_k));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", &(mod_E));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", &(KC));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", &(n_layers));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", &(part_third_centr));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", &(diameter));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", &(Number_crack));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%lf", &(orientation));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", &(parallel));
	fscanf(file, "%s", str);
	fscanf(file, "%s", str);
	fscanf(file, "%u", &(tt));
	fclose(file);
	std::cout <<" s_max = " <<s_max << "; " << " press_k = "<<press_k << " " <<" mod_E = "<< mod_E << " " << " KC = "<<KC << " "
	<< " n_layers = "<<n_layers << " " << "Part = "<<part_third_centr << " " <<"diameter = "<<diameter<< " Number_crack = "<< Number_crack<< " " 
		<< " orientation = "<<orientation<< " "<<" parallel = " <<parallel<< " " <<" tt = "<< tt << " " << std::endl;
	//system("pause");
	/*=====================================*/
	
    //double L_pow = 1;
    //double A_cut = 1.9;

    SetTypeMaterial(0.33); // L_pow = 1;
    SetAngleSprings(0);

    srand((int)time(0));

    area.Periodic() = 1;            // ������������� ��������� �������
    Array<Bond3Dv2> bonds;          // ������ ������ ����� ���������
    Array<AngleBond3D> a_bonds;     // ������ ����� ����� �������
    //int continued = 0;              // ����������� ���������� � �������������� ����� �������� ��������� cnf.dat

    //if (!continued)  { 
		CreateConfig3D(particles, n_layers, area, d_av, bonds, a_bonds); 
	std::cout << "Create config!" << std::endl;
	//}

    Atom3D* as = particles.GetData();
    int n = particles.GetCount();
    int N_B = bonds.GetCount();
    int N_AB = a_bonds.GetCount();
	
	//__________________________________________________________________________________________________________________________________________
	//__________________________________________________________________________________________________________________________________________
	//__________________________________________________________________________________________________________________________________________
	Vector3D center(0.5, 0.5, 0.5);
	Vector3D normal(0, 0, 0);
	normal.y = cos(orientation / 2);
	normal.x = sin(orientation / 2);
	double density = 0.1; //��������� ������������ ������
	//int parallel = 0; // 0 - ������� ��� ���������� ������, 1 - ����������� ����� �� ������
	int intersection = 0; //0 - ��� ����������, 1 - ����
	//CreateFracture_mass(as, area, bonds, density, diameter, parallel, intersection);     // ��������� ������
	std::cout << "normal = " << "(" << normal.x << "," << normal.y << "," << normal.z << ")" << std::endl;
	if ((Number_crack>0.5)&&(Number_crack<2.5)){
		
		CreateFracture_1(as, area, bonds,center, diameter,normal);      // �������� 1 �������
	}else if (Number_crack>3){
		CreateFracture_mass(as, area, bonds, density, diameter, parallel, intersection,orientation);     // ��������� ������
	}
    //CreateFracture_2(as, area, bonds);      // �������� 2� ������
	
	//__________________________________________________________________________________________________________________________________________
	//__________________________________________________________________________________________________________________________________________
	//__________________________________________________________________________________________________________________________________________

    SetB_dt(1.0);                   // ����������� ������

//	SetDifferentModules(area, as, bonds, a_bonds, 1/mod_E, 0);
//	SetDifferentModules(area, as, bonds, a_bonds, mod_E, 1, part_third_centr);
//	SetDifferentModules(area, as, bonds, a_bonds, mod_E, -1, part_third_centr);

//	SetDifferentToughness(area, as, bonds, KC, 1, part_third_centr);
//	SetDifferentToughness(area, as, bonds, KC, -1, part_third_centr);


 /*   if (continued)   {
        LoadConfig("cnf.dat", particles, area, bonds, a_bonds, s0, &press_k_integ);
        as = particles.GetData();
        n = particles.GetCount();
        N_B = bonds.GetCount();
        N_AB = a_bonds.GetCount(); }*/

	/*��������� ����� ������������ �� ���� ����� ������������ ��������� Z = 0*/
	int numConnect = 0;
	int numConnectLiq = 0;
	Bond3Dv2* pB = bonds.GetData();
	Bond3Dv2* pE = pB + bonds.GetCount();
	for (; pB<pE; pB++)           // ���� �� ���� ������
	{
		if (pB->Shifted()) { continue; }
		//if (pB->Gap() < Bond3Dv2::GapM()) { continue; }
		int Firs = pB->First();
		int Sec = pB->Second();
		if (as[Firs].R().z*as[Sec].R().z < 0) {
			numConnect++;
			double rr = area.Height()*0.25*0.5;
			if (pow(as[Firs].R().x, 2) + pow(as[Firs].R().y, 2) < rr*rr)
			{
				numConnectLiq++;
			}
			//if (pB->Gap() < Bond3Dv2::GapM()) { numConnectLiq--; }
		}
		
	}
	std::cout << "ModelSize: " << area.R1().x << "\t" << area.R1().z << "\t" << area.R2().x << "\t" << area.R2().z<<"\n";
	//print("NumCon: %f\t%f\t%f\t%f\n", area.R1().x, area.R1().z, area.R2().x, area.R2().z);
	std::cout << "NumCon: " << numConnect << "\t" << numConnectLiq << "\t" << GetLiq(bonds) << "\n";
	//print("NumCon: %i\t%i\t%i\n", numConnect, numConnectLiq, GetLiq(bonds));
	double d_S = area.Height()*area.Width() / numConnect;
	//print("S: %f\n", area.Height()*area.Width());
	std::cout << "S: " << area.Height()*area.Width() << "\n";
	//system("PAUSE");
	/*----------------------------*/	
	

	/*FILE* f = fopen("output/Config.txt", "w");
	fprintf(f, "%i\t%i\n", 0, bonds.GetCount());
	fclose(f);*/

	int sexit = 100;
	for (s=s0; s<=s_max; s++)
    {

		if (s == 0)
		{
			Get_Init_Coord_N(as, s, area, bonds.GetData(), n, a_bonds.GetData());
		}

        if (outp && !(s%sexit))
        {
            //if (s>s0)   { SaveConfig("cnf.dat", particles, area, bonds, a_bonds, s, &press_k_integ); }
            // ����� ������������ ������
            //sprintf(str, "output/cr%07u.a3r", s);
            //Save_A3R_Gap(str, as, n, bonds, area, 0.5);
			//Save_XYZ_Gap(area,s, as, bonds, n, Volume(bonds, d_S));
			if (s == s0){
				//Save_XYZ_Media(area, s, as, n, mod_E_in, mod_E, mod_E, KC_init,KC, KC, part_third_centr);
				Save_XYZ_Gap(area, s, as, bonds, n, Volume(bonds, d_S));
			}
			if (s%(sexit*4)){
				Save_XYZ_Gap(area,s, as, bonds, n, Volume(bonds, d_S));
			}
			//Stress_XYZ_Gap(as, s, area, bonds.GetData(), N_B, a_bonds.GetData(), n);
			//GetStress_XYZ(area, as, bonds, sigma_x, sigma_y, sigma_z, sigma_xy, sigma_xz, sigma_yz);
			GetStressTensor(as, n, area, bonds.GetData(), N_B, sigma_x, sigma_y, sigma_z, sigma_xy, sigma_xz, sigma_yz);
			//Get_Press(as, s, area, bonds.GetData(), n, N_B, a_bonds.GetData(), d_S, RealP);
			//Get_True_Strain_N(as, s, area, bonds.GetData(), n, a_bonds.GetData(), RealP, diameter * area.Width() / 2);
        }		

        if (outp && !(s%sexit))
        {
           // if (s>s0)   { UpdateParams_E(s_max, press_k, mod_E, KC, n_layers, part_third_centr,diameter); }
            /*print("���: %i\n", s);
            print("����������� ������: %i\n", GetLiq(bonds));
            print("����� �������: %f\n", Volume(bonds,d_S));
            print("�������� ������ ��������: %e\n", press_k);
            print("��������: %e\n\n", press);*/
			std::cout << "Step: " << s
				<< "\nNumber broken bonds: " << GetLiq(bonds)
				<< "\nCrack volume: " << Volume(bonds, d_S)
				<< "\ndQ: " << press_k
				<< "\nPressure: " << press << "\n";

			FILE* f = fopen("output/Crack_Propagation.txt", "a");
			fprintf(f, "%i\t %i\t %e\t\n", s, GetLiq(bonds), Energy_Potencial(bonds, d_S));
            fclose(f);
        }

		/*if (outp && !(s % sexit))
		{
			FILE* fs = fopen("output/Stresses.txt", "a");
			fprintf(fs, "%i\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\t %e\n", s,Energy_Potencial(bonds,d_S), 
			sigma_x, sigma_y, sigma_z, sigma_xy, sigma_xz, sigma_yz, 
			pow(1.000001, s / sexit) - 1, tt*(pow(1-000001, s / sexit) - 1));
			fclose(fs);
		}		*/

		double Vol = Volume(bonds,d_S);
		double press_max = 1.25e-4;

		if ((s % ( v * sexit) == 0))//||(s<20)) 
		{
			//area.ExtendX(as, n, 1-2*0.000001);
			//area.ExtendX(as, n, 1-0.0000001);
		
				// raw compression
				//area.ExtendX(as, n, 1. - 10. * 0.000001);// compression 
				//area.ExtendZ(as, n, 1. + 0.5 * 10 * 0.000001);//  area.ExtendZ(as, n, 1. + 0.33 * 10 * 0.000001);
				//area.ExtendY(as, n, 1. + 0.5 * 10 * 0.000001);//  area.ExtendZ(as, n, 1. + 0.33 * 10 * 0.000001);
			//shear
			std::cout << "Add shear" << std::endl;
			area.ExtendY(as, n, 1.  - 10. * 0.000001);// compression 
			area.ExtendX(as, n, 1.  + 10. * 0.000001);
			
			// shear XY axis

			//area.AddShearYX(as, n , 5. * 0.000001);
			//area.ExtendY(as, n, 1. + 0.5 * 10 * 0.000001);
				//area.ExtendY(as, n, 1-0.0000001*0.5);
				//std::cout << "with shear\n";
			
			//else if (tt < 1.5)
			//{
			//	//area.ExtendY(as, n, 1-2*0.000001*0.25);
			//	area.ExtendX(as, n, 1- 2 * 0.0000001*0.25);
			//}
			
			//area.AddShearYX(as, n, 0.00005);
			//area.AddShearYZ(as, n, 0.0000005);
			//area.AddShearZX(as, n, 0.000001);
			//area.SpinX(as, n, 0.000001);
			//area.SpinY(as, n, 0.000001);
			//area.SpinZ(as, n, 0.000001);
			/*press_k_integ += press_k;
			if (s < 400){
			press = 0.000028;//press_k_integ / (1 + Vol); 
			}
			else {
			press = 0.000028;//press_k_integ / (Vol);
			}*/
			//print("��������: %e\n\n", press);
			//print("����� �������: %f\n", Vol);
		}
        
        //press= min2(press, press_max);
        //press=press_k_integ/Vol;
        //if (press<press_max)    { press_k_integ+=press_k; } //   10^(-3)/100
        //else                    { press=press_max; }

        //if (s<10)   { print("��������: %e\n", press); }
		center = Vector3D(0, 0, 0);
        Step3D(as, n, area, bonds.GetData(), N_B, a_bonds.GetData(), N_AB, press, normal, center, diameter*area.Height());
    }
}
