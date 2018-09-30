#include "atom3Dv2.h"
#include "a3r.h"
#include "tlinear3D.h"
#include "a3r_b.h"
#include "util.h"
//#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <stdio.h>
#include <iostream>
//#include <float.h>
//#include <stdio.h>

///////////////////////////////////////////////////////////////////////////////////////
// определения среднего тензора напряжений по всей области

Tensor GetStressTensor(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, double &(sigma_x), double &(sigma_y), double &(sigma_z), double &(sigma_xy), double &(sigma_xz), double &(sigma_yz))
{
    Tensor tau, tau_s = NULL_T2R;
    int i;
    for (i=0; i<10; i++)
    {
        double sct = 0.2+0.2/3.0*i;
        tau.x = GetStress(as, n, area, bonds, m, 1, sct);//+0.0001);
        tau.y = GetStress(as, n, area, bonds, m, 2, sct);//+0.0002);
        tau.z = GetStress(as, n, area, bonds, m, 3, sct);//+0.0004);
        tau.x += area.ShearYX() * tau.y;
        tau_s += tau;
    }

    tau_s /= 10;

	sigma_x = tau_s.x.x;
	sigma_y = tau_s.y.y;
	sigma_z = tau_s.z.z;
	sigma_xy = tau_s.x.y;
	sigma_xz = tau_s.x.z;
	sigma_yz = tau_s.y.z;

    return tau_s;
	
}

double GetNBonds(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, int axis, double section)
// section д.б. достаточно далеко от 0 и 1, иначе м.б. учтены не все связи (!)
// сдвиг area д.б. =0 (!)
{
    Array<Vector3D> rs(n);
    Atom3D* i;
    Vector3D* j = rs.GetData();
    Atom3D* i_end = as + n;

    for (i = as; i != i_end; i++, j++) //  обратный сдвиг (здесь не учитываем)
    {
        *j = i->R();
        //j->x -= area.ShearYX() * j->y;
    }

    int Axis = axis-1;  // Axis=0,1,2.
    const double* R1 = &area.R1().x;
    const double* R2 = &area.R2().x;
    double w = (1-section)*R1[Axis] + section*R2[Axis];

    int k=0;
    double k_x=0;
    Bond3Dv2* b_i;
    Bond3Dv2* b_end = bonds+m;
    for (b_i=bonds; b_i<b_end; b_i++) // цикл по всем связям
    {
        Vector3D r1(rs[b_i->First()]), r2(rs[b_i->Second()]);
        double* R1_ = &r1.x;
        double* R2_ = &r2.x;
        double  X1 = R1_[Axis];
        double  X2 = R2_[Axis];

        int sort12 = 0;
        if (!(b_i->Shifted() & (3<<(2*Axis))))    // тогда вектор Shift не содержит компоненты вдоль оси Axis
        {                                         // (см. функцию Bond3Dv2::dR(as, area))
            if ((X1 < w) && (X2 >= w)) { sort12 =  1; }
            if ((X2 < w) && (X1 >= w)) { sort12 = -1; }
        }
        if (sort12)
        {
            k++;
            k_x += fabs(X1-X2)/b_i->Length();
        }
    }

    //print("Связей в сечении %i:", axis);
    //print(" %i\n", k);
    //print("Связей на единицу площади: %f\n", k/area.Height()/area.Depth());
    //print("Связей*cos(alpha) на единицу площади: %f\n", k_x/area.Height()/area.Depth());

    //return k/area.Height()/area.Depth();
    return k_x/area.Height()/area.Depth();
    //return k_x/k; == 2/3
}

Vector3D GetStress(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, int axis, double section)
// section д.б. достаточно далеко от 0 и 1, иначе м.б. учтены не все связи (!)
{
    Array<Vector3D> rs(n);
    Atom3D* i;
    Vector3D* j = rs.GetData();
    Atom3D* i_end = as + n;

    for (i = as; i != i_end; i++, j++) //  обратный сдвиг
    {
        *j = i->R();
        j->x -= area.ShearYX() * j->y;
    }

    int Axis = axis-1;  // Axis=0,1,2.
    const double* R1 = &area.R1().x;
    const double* R2 = &area.R2().x;
    double w = (1-section)*R1[Axis] + section*R2[Axis];

    //int k=0;
    Vector3D f = VECT3D_NULL;
    Bond3Dv2* b_i;
    Bond3Dv2* b_end = bonds+m;
    for (b_i=bonds; b_i<b_end; b_i++) // цикл по всем связям
    {
        Vector3D r1(rs[b_i->First()]), r2(rs[b_i->Second()]);
        double* R1_ = &r1.x;
        double* R2_ = &r2.x;
        double  X1 = R1_[Axis];
        double  X2 = R2_[Axis];

        int sort12 = 0;
        if (!(b_i->Shifted() & (3<<(2*Axis))))    // тогда вектор Shift не содержит компоненты вдоль оси Axis
        {                                         // (см. функцию Bond3Dv2::dR(as, area))
            if ((X1 < w) && (X2 >= w)) { sort12 =  1; }
            if ((X2 < w) && (X1 >= w)) { sort12 = -1; }
        }
        if (sort12)
        {
            Vector3D dr(b_i->dR(as, area));
            //double L = dr.Abs();
            double af = b_i->ForceL();

            f.AddProduct(dr, -af*sort12);
            //k++;
        }
    }

    //print("Связей в сечении %i: %i\n", axis, k);
    //nSecBonds++;
    //SecBonds += k;

    f /= area.Volume()/(R2[Axis]-R1[Axis]);
    return f;
}

///////////////////////////////////////////////////////////////////////////////////////

// точка пересечения связи и плоскости трещины
int CrssecBondPlane(Vector3D& crss, Atom3D* as, Bond3Dv2* pB, Vector3D& Center, Vector3D& normal)
{
    Vector3D r1(as[pB->First()].R()), r2(as[pB->Second()].R()); // vector of the particles 1-2 2-1
    r1 -= Center; // r1 = r1 - Center;
    r2 -= Center;
    double  X1 = r1 * normal;
    double  X2 = r2 * normal;
    double  K1 = X2/(X2-X1); //
    double  K2 = X1/(X2-X1);
    crss.Product(K1, r1);
    crss.SubtractProduct(K2, r2);
    return (X1 < 0) != (X2 < 0);
}


// круговая трещина
void CreateFracture(Atom3D* as, Area3d area, Array<Bond3Dv2>& bonds,
                    Vector3D center, double diameter, Vector3D normal)
{
    Vector3D scale = area.R2();
    scale -= area.R1();

    Vector3D Center = scale & center;
    Center += area.R1();

    double Radius2 = sqr(diameter * min3(scale.x, scale.y, scale.z) / 2);
	double a = Radius2;
	double b = Radius2*0.2;
    const char thr=Bond3Dv2::GapM();
    Bond3Dv2* pB = bonds.GetData();
    Bond3Dv2* pE = pB+bonds.GetCount();
    for (; pB<pE; pB++)           // цикл по всем связям
    {
        if (pB->Gap()==0)   { pB->Gap() = 1; }   // все связи разрываемые

        if (!pB->Shifted()) //printf("ERROR!\n"); return;
        {
            Vector3D crss;  // точка пересечения связи и плоскости трещины
            if (CrssecBondPlane(crss, as, pB, Center, normal))
            {
				/*ЭЛЛИПС*/
				/*if ((pow(crss.x, 2) / a + pow(crss.y, 2) / a) < 1) {
					as[pB->First()].N_Gap()++; as[pB->Second()].N_Gap()++;
					pB->Gap() = thr;
				}*/
				/*Квадрат-прямоугольник*/
				/*if ((pow(crss.x, 2) < a) && (pow(crss.y, 2)<b)) {
					as[pB->First()].N_Gap()++; as[pB->Second()].N_Gap()++;
					pB->Gap() = thr;
				}*/
				///*Треугольник*/
				/*if ((crss.y  < crss.x + sqrt(a)) && (crss.y  < -crss.x + sqrt(a)) && (crss.y>-sqrt(a))) {
					as[pB->First()].N_Gap()++; as[pB->Second()].N_Gap()++;
					pB->Gap() = thr;
				}*/

				pB->Gap() = (crss.Sqr() <= Radius2) ? thr : (char)1; 
				
                if (crss.Sqr() <= Radius2)
                {  as[pB->First()].N_Gap()++; as[pB->Second()].N_Gap()++; }
				
            }
        }
    }
}

int Save_A3R_Gap(const char* file_name, Atom3D* as, int n, Array<Bond3Dv2>& bonds, Area3d area, const double r)
{
    // количество разорванных связей при каждой частице
    Array<char> check(n);
    memset(check.GetData(), 0, n);

    // список приграничных частиц
    _List<Atom3D> prt_bond;
    //char thr=6; //  выводятся частицы, у которых не менее 6 разорванных связей
    char thr=1;

    Bond3Dv2* pB = bonds.GetData();
    Bond3Dv2* pE = pB+bonds.GetCount();
    for (; pB<pE; pB++)           // цикл по всем связям
    {
        if (pB->Gap()<Bond3Dv2::GapM()) { continue; } // связь не разорвана
        if (pB->Length() <= 0.9)                      // выводим только короткие связи
        {
            int i=pB->First();
            int j=pB->Second();
            check[i]++;
            check[j]++;
            if (check[i]==thr) { prt_bond.Insert(&as[i]); }
            if (check[j]==thr) { prt_bond.Insert(&as[j]); }
        }
    }

    // добавление окрашенных частиц
    Atom3D* i;
    Atom3D* i_end = as + n;
    for (i=as; i<i_end; i++)                        // цикл по всем частицам
    {
        if (i->Color())     { prt_bond.Insert(i); }
    }

    // добавление угловых частиц области
    Atom3D corner[8];
    for (int i=0; i<8; i++)
    {
        corner[i].R().x = i&1 ? area.R1().x: area.R2().x;
        corner[i].R().y = i&2 ? area.R1().y: area.R2().y;
        corner[i].R().z = i&4 ? area.R1().z: area.R2().z;
        corner[i].Color() = 1;
        prt_bond.Insert(&corner[i]);
    }

    // вывод приграничных частиц
    //return Save_A3R(file_name, prt_bond, r);
    return Save_A3R_B(file_name, prt_bond, r, 0, 1); // используется цвет
}

// количество разорванных связей
int GetLiq(Array<Bond3Dv2>& bonds)
{
    const char thr=Bond3Dv2::GapM();
    int count=0;
    Bond3Dv2* pB = bonds.GetData();
    Bond3Dv2* pE = pB+bonds.GetCount();
    for (; pB<pE; pB++)           // цикл по всем связям
    {
        if (pB->Gap()>=thr) { count++; }
    }
    return count;
}

// объем трещины
//double Volume(Atom3D* as, Area3d area, Array<Bond3Dv2>& bonds)
double Volume(Array<Bond3Dv2>& bonds,double ds)
{
    const char thr=Bond3Dv2::GapM();
    double Vol=0;
    Bond3Dv2* pB = bonds.GetData();
    Bond3Dv2* pE = pB+bonds.GetCount();
    for (; pB<pE; pB++)           // цикл по всем связям
    {
        if (pB->Gap()<thr) { continue; }

        //Vector3D dr(pB->dR(as, area));
        double V0 = pB->dR().Abs()- pB->Length();
        Vol += V0*ds;
        //Vol += V0>0 ? V0 : 0;
    }
    return Vol>0 ? Vol : 0;
}
//Потенциальная энергия

double Energy_Potencial(Array<Bond3Dv2>& bonds,double ds)
{
    const char thr=Bond3Dv2::GapM();
    double E_p=0;
    Bond3Dv2* pB = bonds.GetData();
    Bond3Dv2* pE = pB+bonds.GetCount();
    for (; pB<pE; pB++)           // цикл по всем связям
    {
        if (pB->Gap()<thr) { 
			double dL = pB->dR().Abs()- pB->Length();
			E_p += pB->C_bond()*dL*dL; 
		}

        //Vector3D dr(pB->dR(as, area));
        
        //Vol += V0>0 ? V0 : 0;
    }
    return E_p;
}

//////////////////////////////////////////////////////////////////////////////
// сохранение и загрузка текущей конфигурации частиц для остановки
// и продолжения вычислений с текущего шага по времени

void SaveConfig( const char* file_name, Array<Atom3D>& particles, Area3D& area,
                Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, int& s, double* pV2_max/*=0*/)
{
	FILE *f = fopen(file_name, "wb");
    if (pV2_max==0)  { fwrite("modeling", 8, 1, f); }
    else             { fwrite("model001", 8, 1, f); }

    int n = particles.GetCount();
    void* data = particles.GetData();
    fwrite(&n, 4, 1, f);
    if (n)  { fwrite(data, sizeof(Atom3D), n, f); }

    fwrite(&area, sizeof(Area3D), 1, f);

    n = bonds.GetCount();
    data = bonds.GetData();
    fwrite(&n, 4, 1, f);
    if (n)  { fwrite(data, sizeof(Bond3Dv2), n, f); }

    n = a_bonds.GetCount();
    data = a_bonds.GetData();
    fwrite(&n, 4, 1, f);
    if (n)  { fwrite(data, sizeof(AngleBond3D), n, f); }

    fwrite(&s, 4, 1, f);
    if (pV2_max)  { fwrite(pV2_max, 8, 1, f); }
	fclose(f);
}

void LoadConfig(char* file_name, Array<Atom3D>& particles, Area3D& area,
                Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, int& s, double* pV2_max/*=0*/)
{
	FILE *f = fopen(file_name, "rb");
    char str[10];
    if (pV2_max==0)  { fwrite("modeling", 8, 1, f); }
    else             { fwrite("model001", 8, 1, f); }
    fread(str, 8, 1, f);
    if (pV2_max==0)  {
        if (memcmp(str, "modeling", 8))
        {   print("Ошибка LoadConfig\n");
        return; } }
    else  {
        if (memcmp(str, "model001", 8))
        {   print("Ошибка LoadConfig\n");
        return; } }

    int n;
    void* data;
    fread(&n, 4, 1, f);
    if (n)
    {
        particles.Create(n);
        data = particles.GetData();
        fread(data, sizeof(Atom3D), n, f); }

    fread(&area, sizeof(Area3D), 1, f);

    fread(&n, 4, 1, f);
    if (n)
	{
	bonds.Create(n);
	data = bonds.GetData();
	fread(data, sizeof(Bond3Dv2), n, f); }

	fread(&n, 4, 1, f);
	if (n)
	{
		bonds.Create(n);
		data = a_bonds.GetData();
		fread(data, sizeof(AngleBond3D), n, f);
	}

	fread(&s, 4, 1, f);
	if (pV2_max) { fread(pV2_max, 8, 1, f); }
	fclose(f);
}

void SaveConfig2(char* file_name, Array<Vector3D>& points, int N_equ)
{
	FILE *f = fopen(file_name, "wb");
	fwrite("modeling2", 9, 1, f);

	int n = points.GetCount();
	void* data = points.GetData();
	fwrite(&n, 4, 1, f);
	if (n) { fwrite(data, sizeof(Vector3D), n, f); }

	fwrite(&N_equ, 4, 1, f);
	fclose(f);
}

void LoadConfig2(char* file_name, Array<Vector3D>& points, int& N_equ)
{
	FILE *f = fopen(file_name, "rb");
	char str[10];
	fread(str, 9, 1, f);
	if (memcmp(str, "modeling2", 9))
	{
		print("Ошибка LoadConfig2\n");
		return;
	}

	int n;
	void* data;
	fread(&n, 4, 1, f);
	if (n)
	{
		points.Create(n);
		data = points.GetData();
		fread(data, sizeof(Vector3D), n, f);
	}

	fread(&N_equ, 4, 1, f);
	fclose(f);
}

/***Aspect ratio для трещины****/
double AspectRatio(Atom3D* as, Array<Bond3Dv2>& bonds)
{
	const char thr = Bond3Dv2::GapM();
	double minX = 0, minY = 0, maxX = 0, maxY = 0;
	double Vol = 0;
	Bond3Dv2* pB = bonds.GetData();
	Bond3Dv2* pE = pB + bonds.GetCount();
	for (; pB < pE; pB++)           // цикл по всем связям
	{
		if (pB->Gap() < thr) { continue; }
		int i = pB->First();
		int j = pB->Second();
		double x = 0.5*(as[i].R().x + as[j].R().x);
		double y = 0.5*(as[i].R().y + as[j].R().y);
		if (x > maxX) {
			maxX = x;
		}
		else if (x < minX) {
			minX = x;
		}
		if (y > maxY) {
			maxY = y;
		}
		else if (y < minY) {
			minY = y;
		}
	}
	return (maxX - minX) / (maxY - minY);
}

void SetDifferentModules(Area3D& area, Atom3D* as, Array<Bond3Dv2>& bonds, Array<AngleBond3D>& a_bonds, double factor, int part, double part_third_centr)
{
	double minX = 0, minY = 0, maxX = 0, maxY = 0;
	double Vol = 0;
	Bond3Dv2* pB = bonds.GetData();
	Bond3Dv2* pE = pB + bonds.GetCount();
	double dx = area.Width() / 3;
	int kol_change = 0;
	int kol_not_change = 0;
	double x_area = 0.5*(area.R1() + area.R2()).x;
	for (; pB < pE; pB++)           // цикл по всем связям
	{
		int i = pB->First();
		int j = pB->Second();
		if ((abs(as[i].R().x)< dx*0.49*part_third_centr) && abs(as[j].R().x) < dx*0.49*part_third_centr ){//частицыа в центральном слое
			if (part == 0) {//изменяем центральная слой
				//if ((abs(as[i].R().x - part*dx) < dx*0.49*part_third_centr) && (abs(as[j].R().x - part*dx) < dx*0.49*part_third_centr)) {//обе частицы связи находятся в слое
				pB->Set_C(factor*pB->C_bond());
				pB->Set_Eps(pB->Eps_L() / factor);
				kol_change++;
			}
			else if (part*as[i].R().x > 0) {//изменяем боковые слои
				pB->Set_C(factor*pB->C_bond());
				pB->Set_Eps(pB->Eps_L() / factor);
				kol_change++;
			}
		}
		else {
			kol_not_change++;
		}
	}
	/*проход по угловым пружинкам и увеличение жесткости*/
	AngleBond3D* aB = a_bonds.GetData();
	AngleBond3D* aE = aB + a_bonds.GetCount();
	Vector3D r1, r2, r3; //первый второй третий радиус вектор координат
	for (; aB < aE; aB++) {
		int i = aB->First();
		int j = aB->Second();
		int k = aB->Third();
		r1 = 0.5*(as[bonds[i].First()].R() + as[bonds[i].Second()].R());
		r2 = 0.5*(as[bonds[j].First()].R() + as[bonds[j].Second()].R());
		r3 = 0.5*(as[bonds[k].First()].R() + as[bonds[k].Second()].R());

		if ((abs(r1.x) < dx*0.49*part_third_centr) && abs((r2.x) < dx*0.49*part_third_centr) && abs(r3.x) < dx*0.49*part_third_centr) {
			if (part ==0)
				aB->SetC(factor*aB->GetC());
			else if (part*r1.x > 0) {
				aB->SetC(factor*aB->GetC());
			}
		}
	}
	print("%i\n", kol_change);
	print("%i\n",kol_not_change);
	//print("%i\t%i\n", kol_change, kol_not_change);
}

void SetDifferentToughness(Area3D & area, Atom3D * as, Array<Bond3Dv2>& bonds, double factor, int part, double part_third_centr)
{
	double minX = 0, minY = 0, maxX = 0, maxY = 0;
	double Vol = 0;
	Bond3Dv2* pB = bonds.GetData();
	Bond3Dv2* pE = pB + bonds.GetCount();
	double dx = area.Width() / 3;
	int kol_change = 0;
	int kol_not_change = 0;
	for (; pB < pE; pB++)           // цикл по всем связям
	{
		int i = pB->First();
		int j = pB->Second();
		if ((abs(as[i].R().x)< dx*0.49*part_third_centr) && abs(as[j].R().x) < dx*0.49*part_third_centr) {//связь в центральном слое
			if (part == 0) {//изменяем центральная слой
				pB->Set_Eps(factor*pB->Eps_L());
				kol_change++;
			}
			else if (part*as[i].R().x > 0) {//изменяем боковые слои
				pB->Set_Eps(factor*pB->Eps_L());
				kol_change++;
			}
		}
		else {
			kol_not_change++;
		}
	}
}
/*Функции сохранения в файл для отображения*/

void Save_XYZ_Media(Area3D & area, int n_step, Atom3D * as, int n, double E_centr,double E_left,double E_right, double KC_centr,double KC_left, double KC_right,double part_third)
{
	const char thr = Bond3Dv2::GapM();
	int count = 0;
	char str[100];
	sprintf(str, "output/Media%07u.xyz", n_step);
	int number = 0;
	std::vector<int> numPartGap;
	double dx = area.Width() / 3;
	for (int i = 0; i < n; i++)
	{
		if (as[i].N_Gap() >= thr) { continue; }
		number++;
		numPartGap.push_back(i);
	}
	FILE* f = fopen(str, "w");
	fprintf(f, "%i\n", number);
	fprintf(f, "%i\n", number);
	int k = 1;
	for (int it =0;it<numPartGap.size();it++)
	{
		int i = numPartGap.at(it);
		if (abs(as[i].R().x) <= dx*0.5*part_third)//
		//if ((abs(as[i].R().x - part*dx) < dx*0.49*part_third_centr){
		{
			fprintf(f, "%i %f %f %f %f %f %\n", k, as[i].R().x, as[i].R().y, as[i].R().z, E_centr, KC_left);
		}
		else if (as[i].R().x >0)
		{
			fprintf(f, "%i %f %f %f %f %f %\n", k, as[i].R().x, as[i].R().y, as[i].R().z, E_right, KC_right);
		}
		else {
			fprintf(f, "%i %f %f %f %f %f %\n", k, as[i].R().x, as[i].R().y, as[i].R().z, E_left, KC_left);
		}

	}
	fclose(f);
}


void Stress_XYZ_Gap(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds, int num)
{
	const char thr = Bond3Dv2::GapM();
	int number = 0;
	Bond3Dv2* b_i;
	Bond3Dv2* b_end = bonds + m;
	

	for (b_i = bonds; b_i < b_end; b_i++)
	{
		if ((as[b_i->First()].N_Gap() < thr) && 
			(as[b_i->Second()].N_Gap() < thr) && 
			(b_i->Gap() != thr))
			number++;
	}

	std::vector<double> stress(num);

	char str[100];
	sprintf(str, "output/stress/stress%07u.xyz", n_step);


	for (b_i = bonds; b_i < b_end; b_i++)               // цикл по всем связям
	{
		double L = b_i->dR().Abs();         //  значение dR сохраняется в переменной (*b_i)
		double L0 = b_i->Length();
													// условие разрыва связи
		int i = 0;
		int j = 0;
		if ((as[b_i->First()].N_Gap() < thr) &&
			(as[b_i->Second()].N_Gap() < thr) &&
			(b_i->Gap() != thr))
		{
			i = b_i->First();
			j = b_i->Second();
			stress[i] += L * b_i->ForceL() / (2 * area.Volume() / num); //* b_i->dR() / b_i->dR().Abs();
			stress[j] += L * b_i->ForceL() / (2 * area.Volume() / num); //* b_i->dR() / b_i->dR().Abs();
		}

	}

	int j = 0;

	for (int i = 0; i < stress.size(); i++)
		if (abs(stress[i]) != 0)
			j += 1;

	FILE* f = fopen(str, "w");
	fprintf(f, "%i\n", j);
	fprintf(f, "Volume  \t%f\n", 0);
	int k = 1;

	for (int i = 0; i < stress.size(); i++)
		if (abs(stress[i]) != 0)
			fprintf(f, "%i %f %f %f %f %\n", k, as[i].R().x, as[i].R().y, as[i].R().z, stress[i]);

	fclose(f);
}

void Strain_XYZ_Gap(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds)
{
	const char thr = Bond3Dv2::GapM();
	int number = 0;
	Bond3Dv2* b_i;
	Bond3Dv2* b_end = bonds + m;

	for (b_i = bonds; b_i < b_end; b_i++)
	{
		if (b_i->Gap() == thr)
			number++;
	}
	
	char str[100];
	char stre[100];
	sprintf(str, "output/strain/strain%07u.xyz", n_step);
	sprintf(stre, "output/strain/Раскрытие%07u.txt", n_step);

	FILE* f = fopen(str, "w");
	FILE* file = fopen(stre, "w");

	fprintf(f, "%i\n", number);
	fprintf(f, "Volume  \t%f\n", 0);


	
	for (b_i = bonds; b_i < b_end; b_i++)               // цикл по всем связям
	{
		double L = b_i->dR(as, area).Abs();         //  значение dR сохраняется в переменной (*b_i)
		double L0 = b_i->Length();
		double _eps = L0 / L - 1;                   //  деформация длины связи (>0 при сжатии)
													// условие разрыва связи
		int i = 0;
		int j = 0;

		if (b_i->Gap() == thr) //- b_i->Eps_L()
		{
			i = b_i->First();
			j = b_i->Second();
			//fprintf(f, "%i %f %f %f %f %i %i %i %\n", k, as[i].R().x, as[i].R().y, as[i].R().z, L - L0, 1, 0, 0);
			//fprintf(file, "%f %f %f %\n", as[i].R().x, as[i].R().y, as[i].R().z);
			//if (as[i].R().x * as[j].R().x < 0)
				fprintf(f, "%i %f %f %f %f %i %i %i %\n", 1, as[i].R().x, as[i].R().y, as[i].R().z, L - L0, 1, 0, 0);
		}
	}
	
	fclose(f);
	fclose(file);
	
}

void Save_XYZ_Gap(Area3D & area, int n_step, Atom3D* as, Array<Bond3Dv2>& bonds, int n, double dV)
{
	const char thr = Bond3Dv2::GapM();
	int number = 0;
	char str[100];
	sprintf(str, "output/crack%07u.xyz", n_step);
	int count = 0;
	std::vector<int> numPartGap;
	
	for (int i = 0; i < n; i++)
	{
		if (as[i].N_Gap() < thr) { continue; }
		number++;
		numPartGap.push_back(i);
	}
	Bond3Dv2* pB = bonds.GetData();
	Bond3Dv2* pE = pB + bonds.GetCount();
	//double d_V = 0;
	//for (; pB<pE; pB++)           // цикл по всем связям
	//{
	//	if (pB->Shifted()) { continue; }
	//	//if (pB->Gap() < Bond3Dv2::GapM()) { continue; }
	//	if (pB->Gap() == thr)
	//	{
	//		int i = pB->First();
	//		int j = pB->Second();
	//		double L = pB->dR(as, area).Abs();         //  значение dR сохраняется в переменной (*b_i)
	//		d_V += abs(L - pB->Length())*d_s;
	//	}
	//}
	//	
	FILE* f = fopen(str, "w");
	fprintf(f, "%i\n", number);
	fprintf(f, "Volume  \t%f\n", dV);
	int k = 1;
	for (int it = 0; it<numPartGap.size(); it++)
	{
		int i = numPartGap.at(it);
		fprintf(f, "%i %f %f %f %f %\n", k, as[i].R().x, as[i].R().y, as[i].R().z, 1.0);
	}
	fclose(f);
}

void GetStress_XYZ(Area3d area, Atom3D* as, Array<Bond3Dv2>& bonds, double &(sigma_x), double &(sigma_y), double &(sigma_z), double &(sigma_xy), double &(sigma_xz), double &(sigma_yz))
{
	const char thr = Bond3Dv2::GapM();
	double sum_sigma_x = 0;
	double sum_sigma_y = 0;
	double sum_sigma_z = 0;
	Bond3Dv2* pB = bonds.GetData();
	Bond3Dv2* pE = pB + bonds.GetCount();
	Vector3D tau_x(0, 0, 0);
	Vector3D tau_y(0, 0, 0);
	Vector3D tau_z(0, 0, 0);
	double d_av = 1;
	for (; pB < pE; pB++)           // цикл по всем связям
	{
		int i = pB->First();
		int j = pB->Second();
		double forse = pB->ForceL();
		Vector3D e1 = pB->dR(as,area) / pB->dR(as,area).Abs();
		if ((as[i].R().z*as[j].R().z<0) && (abs(as[i].R().z)<area.Depth()*0.2))
		{
			if (as[i].R().z > 0)
			{
				tau_z += forse*e1;
			}
			else
			{
				tau_z -= forse*e1;
			}
		}
		else if ((as[i].R().x*as[j].R().x < 0) && (abs(as[i].R().x)<area.Height()*0.2))
		{
			if (as[i].R().x > 0)
			{
				tau_x += forse*e1;
			}
			else
			{
				tau_x -= forse*e1;
			}
		}
		else if ((as[i].R().y*as[j].R().y < 0) && (abs(as[i].R().y)<area.Width()*0.2))
		{
			if (as[i].R().y > 0)
			{
				tau_y += forse*e1;
			}
			else
			{
				tau_y -= forse*e1;
			}
		}
	}

	double VV = (area.Height())* (area.Width()) * 1;
	sigma_x = tau_x.x / VV;
	sigma_y = tau_y.y / VV;
	sigma_z = tau_z.z / VV;
	sigma_xy = tau_y.x / VV;
	sigma_xz = tau_z.x / VV;
	sigma_yz = tau_y.z / VV;
	//std::cout << "TAU X: " << tau_x.x << " " << tau_x.y << " " << tau_x.z << "\n";
	//std::cout << "TAU Y: " << tau_y.x << " " << tau_y.y << " " << tau_y.z << "\n";
	//std::cout << "TAU Z: " << tau_z.x << " " << tau_z.y << " " << tau_z.z << "\n";

	
}

void Get_True_Strain(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds)
{
	const char thr = Bond3Dv2::GapM();
	int number = 0;
	Bond3Dv2* b_i;
	Bond3Dv2* b_end = bonds + m;

	for (b_i = bonds; b_i < b_end; b_i++)
	{
		if (b_i->Gap() == thr)
			number++;
	}

	char str[100];
	//char stre[100];
	sprintf(str, "output/True_Strain/TrueStrain%07u.xyz", n_step);
	//sprintf(stre, "output/Maximum.txt", n_step);

	FILE* f = fopen(str, "w");
	FILE* fdop = fopen("output/Maximum.txt", "a");

	fprintf(f, "%i\n", number);
	fprintf(f, "Volume  \t%f\n", 0);

	double zmax = 0;
	double xmax = 0;


	for (b_i = bonds; b_i < b_end; b_i++)               // цикл по всем связям
	{
		double L = b_i->dR(as, area).Abs();         //  значение dR сохраняется в переменной (*b_i)
		double L0 = b_i->Length();
		double _eps = L0 / L - 1;                   //  деформация длины связи (>0 при сжатии)
													// условие разрыва связи
		int i = 0;
		int j = 0;

		if (b_i->Gap() == thr) //- b_i->Eps_L()
		{
			i = b_i->First();
			j = b_i->Second();
			
			fprintf(f, "%i %e %e %e %e\n", 1, as[i].R().x, as[i].R().y, as[i].R().z, abs(as[i].R().z - b_i->InitZ1()));
			
			if (abs(as[i].R().z - b_i->InitZ1()) >= zmax)
				zmax = abs(as[i].R().z - b_i->InitZ1());

			if (abs(as[i].R().x) >= xmax)
				xmax = abs(as[i].R().x);
		}
	}

	fprintf(fdop, "%e %e %\n", xmax, zmax);

	fclose(f);
	fclose(fdop);
}

void Get_Init_Coord(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds)
{
	const char thr = Bond3Dv2::GapM();
	Bond3Dv2* b_i;
	Bond3Dv2* b_end = bonds + m;

	for (b_i = bonds; b_i < b_end; b_i++)               // цикл по всем связям
	{
			b_i->Set_InitZ1(as[b_i->First()].R().z);
			b_i->Set_InitZ2(as[b_i->Second()].R().z);
	}
}

void Get_True_Strain_N(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int n, AngleBond3D* a_bonds, double P, double R)
{
	const char thr = Bond3Dv2::GapM();

	int number = 0;
	//double FA = 0;
	//double FN = 0;

	for (int i = 0; i < n; i++)
	{
		if ((as[i].N_Gap() < thr) || (as[i].R().z < 0)) { continue; }
		number++;
	}

	char str[100];
	//char stre[100];
	sprintf(str, "output/True_Strain/TrueStrain%07u.xyz", n_step);
	//sprintf(stre, "output/Maximum.txt", n_step);

	FILE* f = fopen(str, "w");
	//FILE* fdop = fopen("output/Maximum.txt", "a");

	fprintf(f, "%i\n", number);
	fprintf(f, "Volume  \t%f\n", 0);

	//double zmax = 0;
	//double xmax = 0;

	for (int i = 0; i < n; i++)
	{
		if (as[i].N_Gap() < thr) { continue; }

		if (as[i].R().z >= 0)
			fprintf(f, "%i\t %e\t %e\t %e\t %e\n", 1, as[i].R().x, as[i].R().y, as[i].R().z, abs(as[i].R().z - as[i].Init_Z()));

		/*FN += abs(as[i].R().z - as[i].Init_Z()) / number;

		if (1 - pow((as[i].R().x / R), 2) - pow((as[i].R().y / R), 2) > 0)
			FA += 8 * P * R * (1 - 0.256031 * 0.256031) * sqrt(1 - pow((as[i].R().x / R), 2) - pow((as[i].R().y / R), 2)) / 6.41 / PI / 2 / number;*/
		/*if (abs(as[i].R().z - as[i].Init_Z()) >= zmax)
			zmax = abs(as[i].R().z - as[i].Init_Z());

		if (abs(as[i].R().x) >= xmax)
			xmax = abs(as[i].R().x);*/
		
	}	

	//fprintf(fdop, "%i\t %e\t %e\t %\n", n_step, FN, FA);

	fclose(f);
	//fclose(fdop);
}

void Get_Init_Coord_N(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int n, AngleBond3D* a_bonds)
{
	const char thr = Bond3Dv2::GapM();

	for (int i = 0; i < n; i++)
	{
		as[i].Set_Init_R(as[i].R()); // initial particle position
		if (as[i].N_Gap() < thr) { continue; }
		as[i].Set_Init_Z(as[i].R().z);
		
	}

}

void Get_Press(Atom3D* as, int n_step, Area3d area, Bond3Dv2* bonds, int n, int m, AngleBond3D* a_bonds, double d_S, double &(RealP))
{
	const char thr = Bond3Dv2::GapM();
	int number_bond = 0;
	int number_part = 0;
	double F_sum = 0;
	double F_sum_norm = 0;
	double F_sum_L = 0;
	double F_sum_norm_L = 0;
	Bond3Dv2* b_i;
	Bond3Dv2* b_end = bonds + m;

	for (b_i = bonds; b_i < b_end; b_i++)
	{
		if (b_i->Gap() == thr)
			number_bond++;
	}

	for (int i = 0; i < n; i++)
	{
		if (as[i].N_Gap() < thr) { continue; }
			number_part++;
	}

	for (b_i = bonds; b_i < b_end; b_i++)               // цикл по всем связям
	{
		if (b_i->Gap() == thr)
		{
			Vector3D z;
			z.Set(0, 0, 1);
			F_sum += b_i->ForceL();
			F_sum_norm += b_i->ForceL() * abs((z * (as[b_i->First()].R() - as[b_i->Second()].R())) / (z.Abs() * (as[b_i->First()].R() - as[b_i->Second()].R()).Abs()));

			F_sum_L += b_i->ForceL() * b_i->dR().Abs();
			F_sum_norm_L += b_i->dR().Abs() * b_i->ForceL() * abs((z * (as[b_i->First()].R() - as[b_i->Second()].R())) / (z.Abs() * (as[b_i->First()].R() - as[b_i->Second()].R()).Abs()));
		}
	}

	//char str[100];
	//sprintf(str, "output/Real_Press.txt", n_step);

	FILE* f = fopen("output/Real_Press.txt", "a");
	fprintf(f, "%e\t %e\t %e\t %e\n", F_sum / number_bond, F_sum_norm / number_bond, F_sum / (d_S * number_bond), F_sum_norm / (d_S * number_bond));
	fclose(f);

	//char stre[100];
	//sprintf(stre, "output/Real_Press_L.txt", n_step);

	FILE* fe = fopen("output/Real_Press_L.txt", "a");
	fprintf(fe, "%e\t %e\t %e\t %e\n", F_sum_L / number_bond, F_sum_norm_L / number_bond, F_sum_L / (d_S * number_bond), F_sum_norm_L / (d_S * number_bond));
	fclose(fe);

	RealP = F_sum_norm_L / (d_S * number_bond);
}