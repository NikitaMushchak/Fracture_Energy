//-----------------------------------------------------------------------------------

//	Класс 3D векторов

//	from 'vect3d.h' /Anton M. Krivtsov/  2001

//	В. Цаплин

//	2008

//  Изменен 19.06.2009
//  Изменен 16.09.2011

//-----------------------------------------------------------------------------------

#include "vector3d.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>

//-----------------------------------------------------------------------------------

//static const double PI4 = PI * 4;

//-----------------------------------------------------------------------------------

const double* pAMax2(const double& x1, const double& x2)
{
	return (fabs(x1) > fabs(x2)) ? &x1 : &x2;
}

//-----------------------------------------------------------------------------------
//	RAND_MAX = 32767 = 2^15 - 1
//-----------------------------------------------------------------------------------

double Vector3D::Abs() const	{ return sqrt(Sqr()); }

//-----------------------------------------------------------------------------------

void Vector3D::SetRand(double x_min, double y_min, double z_min, double x_max, double y_max, double z_max) 
{
	x = rand(x_min, x_max);
	y = rand(y_min, y_max);
	z = rand(z_min, z_max);
}

//-----------------------------------------------------------------------------------

void Vector3D::SetRand(double r_max)
{
	double rr_max = r_max * r_max;
	do SetRand(-r_max, -r_max, -r_max, r_max, r_max, r_max);
	while (Sqr() > rr_max);
}

//-----------------------------------------------------------------------------------

double Vector3D::GetNorm1() const	{ return max3( fabs(x), fabs(y), fabs(z) ); }

//-----------------------------------------------------------------------------------

void Vector3D::SetTurn(double a, Vect3d axis)
{
	double b = a / PI4;	b -= floor(b);

	double t;

	if (fabs(b - 0.5) > 0.25)
	{
		t = tan(a / 4);
	}
	else
	{
		t = tan((a - PI2) / 4);
	}

	x = axis.x * t; y = axis.y * t; z = axis.z * t;
}

//-----------------------------------------------------------------------------------

double Vector3D::GetTurnAngle()	const
{
	return 4 * atan(Abs());
}

//-----------------------------------------------------------------------------------
// добавление дополнительного малого поворота
// устанавливает безопасный вектор поворота (Abs() <= sqrt(2))

void Vector3D::AddSmallTurn(Vect3d dt)
{
	/*Vector3D _2dt;

	double a = Sqr(), b = dt * (*this);

	_2dt.Product(dt, (1 - a) / 2);
	_2dt.AddVectorProduct(dt, *this);
	_2dt.AddProduct(*this, b);

	x += _scalb(_2dt.x, -1);
	y += _scalb(_2dt.y, -1);
	z += _scalb(_2dt.z, -1);

	a += (1 + a) / 2 * b;

	if (a > 2)   
	{
		x /= -a; y /= -a; z /= -a;
	}*/
}

//-----------------------------------------------------------------------------------
// добавление дополнительного поворота
// устанавливает нормализованный вектор поворота (Abs() <= 1)

void Vector3D::AddTurn(Vect3d t)
{
	/*Vector3D t_(*this);

	double t12 = t_.Sqr();
	double t22 = t.Sqr();
	double p = t * t_;
	
	p = _scalb(p, 1);

	double A = 1 + t12 * t22 - p;
	double B = t12 + t22 + p;

	if (A < B)   A = -B;

	Product   (t_, (1 - t22) / A);
	AddProduct(t,  (1 - t12) / A);

	p = 2 / A; t_.x *= p; t_.y *= p; t_.z *= p;

	AddVectorProduct(t, t_);*/
}

//-----------------------------------------------------------------------------------

void Vector3D::NormalizeTurn()
{
	double a = Sqr();

	if (a > 1)
	{
		x /= -a; y /= -a; z /= -a;
	}
}

//-----------------------------------------------------------------------------------
// устанавливает повернутый вектор v поворотом t

void Vector3D::Rotate(Vect3d t, Vect3d v)
{
	double t2 = t.Sqr(), _1t2 = t2 + 1, _1_t2 = 1 - t2, _1t22 = _1t2 * _1t2;
	//double A = (_1_t2 * _1_t2 - 4 * t2) / _1t22;
	double A = (1 + t2 * (t2 - 6)) / _1t22;
	double B = 8 * (t * v) / _1t22;
	double C = 4 * _1_t2 / _1t22;

	Vector3D tv; tv.VectorProduct(t, v);

	Product(A, v);
	AddProduct(B, t);
	AddProduct(C, tv);
}

//-----------------------------------------------------------------------------------
