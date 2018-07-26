#include "atom3Dv2.h"
//#include "util.h"
#define _USE_MATH_DEFINES
#include <math.h>
//#include <float.h>
//#include <stdio.h>

static double d_min_d=1;
static double T0 = M_PI*d_min_d;
static double dt = 0.02*T0;
static const double C_bond = 1;
static const double om = sqrt(C_bond);
//static double b_dt = 10 * om * dt; // om * dt = критическое трение
static double b_dt = om * dt; // om * dt = критическое трение

void Set_d_min_d_step(double d_min_d, double L_pow)
{
    ::d_min_d=d_min_d;
    T0 = M_PI*pow(d_min_d, L_pow/2); // на основе ГЦК
    dt = 0.02*T0;
}

void SetB_dt(double k) { b_dt=k * om * dt; }  // устанавливаем коэффициент трения


///////////////////////////////////////////////////////////////////////////////////////

//  вычисление среднего перемещения частиц за один шаг

double GetLastAvDisplacement(Atom3D* as, int n)
{
    double D = 0.0;
    Atom3D* i;
    Atom3D* i_end = as + n;

    for (i = as; i != i_end; i++)
    {
        D += i->V().Abs();
        //D += i->V().Sqr();
    }

    return D/n*dt;
    //return sqrt(D/n)*dt;
}

//  вычисление среднего перемещения частиц относительно начальных положений

double GetFullAvDisplacement(Atom3D* as, int n, Vector3D* points0)
{
    double sum=0;
    Atom3D* i;
    Vector3D* j = points0;
    Atom3D* i_end = as + n;

    for (i = as; i != i_end; i++, j++)
    {
        Vector3D D;
        D.Difference(i->R(), *j);
        sum += D.Abs();
    }

    return sum/n;
}

double GetFullAvDisplacement(Atom3D* as, int n, Atom3D* points0)
{
    double sum=0;
    Atom3D* i;
    int j = 0;
    Atom3D* i_end = as + n;

    for (i = as; i != i_end; i++, j++)
    {
        Vector3D D;
        D.Difference(i->R(), points0[j].R());
        sum += D.Abs();
    }

    return sum/n;
}

double GetFullAvDisplacementGap(Atom3D* as, int n, Vector3D* points0)
{
    double sum=0;
    Atom3D* i;
    Vector3D* j = points0;
    Atom3D* i_end = as + n;
    int m=0;

    for (i = as; i != i_end; i++, j++)
    if (i->N_Gap())
    {
        Vector3D D;
        D.Difference(i->R(), *j);
        sum += D.Abs();
        m++;
    }

    return sum/m;
}

///////////////////////////////////////////////////////////////////////////////////////
// система с угловыми пружинками

/*
void Step3D_angle(Bond3Dv2* bonds, AngleBond3D* a_bonds, int p)
{
    const char thr=Bond3Dv2::GapM();
    AngleBond3D* ab_i;
    AngleBond3D* ab_end = a_bonds+p;
    for (ab_i=a_bonds; ab_i<ab_end; ab_i++) // цикл по всем углам
    {
        Bond3Dv2& b1 = bonds[ab_i->First()];
        Bond3Dv2& b2 = bonds[ab_i->Second()];
        if (b1.Gap()==thr || b2.Gap()==thr) { continue; } // одна из связей, прилежащих к углу, разорвана
        int i=b1.First();
        int j=b1.Second();
        Vector3D& R1 = b1.dR();
        Vector3D& R2 = b2.dR();
        double R1_2 = R1.Sqr();
        double R2_2 = R2.Sqr();
        double R1R2 = R1*R2;
        double vR = sqrt(R1_2*R2_2 - R1R2*R1R2); // vR=|R1 X R2|
        if (!(vR>0))  continue;  // если vR==0 или неопределенность, то силы==0
        int sqn1 = (b2.First()==i || b2.Second()==i) ? 1 : -1;
        int sgn2 = (b2.First()==i || b2.First()==j) ? 1 : -1;
        int sgn = sgn1*sgn2;
        double dfi = acos(R1R2*sgn/sqrt(R1_2*R2_2)) - ab_i->Angle();

        Vector3D dV1, dV2;
        dV1.Product(R1R2/R1_2, R1);
        dV2.Product(R1R2/R2_2, R2);
        dV1 -= R2;
        dV2 -= R1;
        dV1 *= c_fi*dfi/vR*dt*sgn;//?
        dV2 *= c_fi*dfi/vR*dt*sgn;//?

        as[i] += dV1;
        as[j] -= dV1;
        as[b2.First()]  += dV2;
        as[b2.Second()] -= dV2;
    }
}
*/

static double c_fi=0;                       // жесткость угловой пружинки

void   SetAngleSprings(double C_fi, double a_cut1/*=0*/, double a_cut2/*=0*/)
{
    c_fi = C_fi;
    SetAcut12(a_cut1, a_cut2);
}
double GetC_fi()            { return c_fi; }

void Step3D_angle(Bond3Dv2* bonds, AngleBond3D* a_bonds, int p)
{
    const char thr=Bond3Dv2::GapM();
    AngleBond3D* ab_i;
    AngleBond3D* ab_end = a_bonds+p;
    for (ab_i=a_bonds; ab_i<ab_end; ab_i++)             // цикл по всем углам
    {
        Bond3Dv2& b1 = bonds[ab_i->First()];
        Bond3Dv2& b2 = bonds[ab_i->Second()];
        Bond3Dv2& b3 = bonds[ab_i->Third()];            // сторона треугольника, противолежащая углу
        if (b1.Gap()==thr || b2.Gap()==thr) { continue; } // одна из связей, прилежащих к углу, разорвана
        Vector3D& R1 = b1.dR();
        Vector3D& R2 = b2.dR();
        double R1_2 = R1.Sqr();
        double R2_2 = R2.Sqr();
        double R1R2 =
            (b2.First()==b1.First() ||
             b2.Second()==b1.Second()) ?                // если направление R1 или R2 д.б. изменено,
            (R1*R2): -(R1*R2);                          // меняем знак R1*R2
        double vR = sqrt(R1_2*R2_2 - R1R2*R1R2);        // ==|R1|*|R2|*sin(fi)
        if (!(vR>0)) { continue; }                      // если vR==0 или неопределенность, то силы==0

        double k_dfi = atan2(vR, R1R2) - ab_i->Angle(); // изменение угла (delta_fi)
      //k_dfi *= c_fi/vR;
        k_dfi *= ab_i->GetC()/vR;

        b1.ForceL() += float(k_dfi * (1 - R1R2/R1_2));  //  добавляем силы
        b2.ForceL() += float(k_dfi * (1 - R1R2/R2_2));  //  от угловых
        b3.ForceL() -= float(k_dfi);                    //  пружинок
    }
}

// press - величина, пропорциональная давлению внутри трещин (сила взаим. в разорванных связях)
// Давление задается правее координаты X, равной (*pFacet)

void Step3D(Atom3D* as, int n, Area3d area, Bond3Dv2* bonds, int m, AngleBond3D* a_bonds, int p,
            double press/*=0*/, double* pFacet/*=0*/, _List<Atom3D>* pListSlidingSealingY/*=0*/)
{
    const char thr=Bond3Dv2::GapM();
    Bond3Dv2* b_i;
    Bond3Dv2* b_end = bonds+m;

    double Facet;                                   // грань, правее которой задается давление
    if (pFacet)  { Facet=area.X1() + (*pFacet)*area.Width(); }

    for (b_i=bonds; b_i<b_end; b_i++)               // цикл по всем связям
    {
        double L = b_i->dR(as, area).Abs();         //  значение dR сохраняется в переменной (*b_i)
        double L0 = b_i->Length();
        double _eps = L0 / L - 1;                   //  деформация длины связи (>0 при сжатии)
        double f_L = b_i->C_bond() * _eps;          // C_bond = C1;

        // условие разрыва связи
		if (b_i->Gap() && b_i->Gap() < thr)           // разрываемая, но не разорванная связь
		{
			//if (f_L<-0.0002)                      // условие разрыва
			//if (f_L*L0<-0.0003)
			if (_eps < -b_i->Eps_L())
			{
				b_i->Gap()++;   // разрываем связь
				/*Характеризуем частицы с порванными связями*/
				int i = b_i->First();
				int j = b_i->Second();
				as[i].N_Gap()++;
				as[j].N_Gap()++;

			}
        //else                 { b_i->Gap()=1; }  // не разорванная связь
        }

        // сила отталкивания / L (давление или продольные пружинки)
        if (b_i->Gap()==thr)
        {
			//kontact here
            //b_i->ForceL() = (!pFacet || as[b_i->First()].R().x>Facet) ? float(press / L): 0;
			/*if (L  < L0)
				b_i->ForceL() = float(2*f_L);*/
			b_i->ForceL() = ( L < L0 ) ? float(10 * f_L) : 0;
        }
        else        
		{ 
			b_i->ForceL() = float(f_L);
			/*
			else { b_i->ForceL() = 0; }*/
		}
    }

    if (p)  { 
		Step3D_angle(bonds, a_bonds, p); 
	}    // силы от угловых пружинок

    for (b_i=bonds; b_i<b_end; b_i++)               // цикл по всем связям
    {
        Vector3D dV;
        dV.Product(b_i->dR(), b_i->ForceL() * dt);
        as[b_i->First()]  -= dV;                    // интегрирование
        as[b_i->Second()] += dV;                    // ускорения
    }

    if (pListSlidingSealingY)                       // имеется скользящая заделка на грани, ортог. оси OY
    {
        BOOL more = pListSlidingSealingY->Restore();
        while(more)
        {
            //Atom3D* i = pListSlidingSealingY->Iterate(more);
            //i->V().y  = 0;
            pListSlidingSealingY->Iterate(more)->V().y  = 0;
        }
    }

    Atom3D* i;
    Atom3D* i_end = as + n;
    for (i=as; i<i_end; i++)                        // цикл по всем частицам
    {
        i->V().SubtractProduct(0.225 * b_dt, i->V());       // вязкое трение
        i->R().AddProduct(i->V(), dt);              // интегрирование скорости
    }
}
