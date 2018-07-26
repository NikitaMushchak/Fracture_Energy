//_________________________________________________________________________________________

//	 ласс 3D векторов

//	from 'vect3d.h' /Anton M. Krivtsov/  2001

//	¬. ÷аплин

//	»зменен 2008

//  »зменен 19.06.2009
//  »зменен 14.09.2011
//  »зменен 01.12.2016

//_________________________________________________________________________________________

#ifndef ___Vector3D_h___
#define ___Vector3D_h___

//_________________________________________________________________________________________

#include <stdlib.h>

static const double PI = 3.1415926535897932;
static const double PI2 = PI * 2;
static const double PI4 = PI * 4;

//_________________________________________________________________________________________

#ifndef sqr
#define sqr square
template <class T> inline T sqr(const T& x) { return x * x; }
#endif

template <class T> inline T cube(const T& x) { return x * x * x; }

//_________________________________________________________________________________________

template <class T> inline T max2(const T& x1, const T& x2) { return (x1 > x2) ? x1 : x2; }
template <class T> inline T min2(const T& x1, const T& x2) { return (x1 < x2) ? x1 : x2; }

//_________________________________________________________________________________________

template <class T> inline const T* pMax2(const T& x1, const T& x2)
{ return (x1 > x2) ? &x1 : &x2; }

const double* pAMax2(const double& x1, const double& x2);

//_________________________________________________________________________________________

template <class T> inline T max3(const T& x1, const T& x2, const T& x3) 
{ return (x1 > x2) ? max2(x1, x3) : max2(x2, x3); }

template <class T> inline T min3(const T& x1, const T& x2, const T& x3) 
{ return (x1 < x2) ? min2(x1, x3) : min2(x2, x3); }

//_________________________________________________________________________________________

template <class T> inline const T* pMax3(const T& x1, const T& x2, const T& x3)
{ return (x1 > x2) ? pMax2(x1, x3) : pMax2(x2, x3); }

//_________________________________________________________________________________________

template <class T> inline T max4(const T& x1, const T& x2, const T& x3, const T& x4) 
{ return max2( max2(x1, x2), max2(x3, x4) ); }

template <class T> inline T min4(const T& x1, const T& x2, const T& x3, const T& x4) 
{ return min2( min2(x1, x2), min2(x3, x4) ); }

//_________________________________________________________________________________________

inline const double* pAMax4(const double& x1, const double& x2,
							const double& x3, const double& x4) 
{ return pAMax2( *pAMax2(x1, x2), *pAMax2(x3, x4) ); }

//_________________________________________________________________________________________

inline double rand(double min, double max)	//	случайное число от min до max
{ 
	return min + (max - min) * rand() / RAND_MAX; 
}

//_________________________________________________________________________________________

#define BOOL	int

#define Vect3d const Vector3D&

//_________________________________________________________________________________________

class Vector3D
{
public:
	Vector3D()									{}				
	Vector3D(double x, double y, double z)		{ Set(x, y, z); } 
	Vector3D(Vect3d  r)							{ x = r.x; y = r.y; z = r.z; }

	void Set(double x, double y, double z)		{ this->x = x; this->y = y; this->z = z; }

	void SetRand(double x_min, double y_min, double z_min, double x_max, double y_max, double z_max);
	void SetRand(double r_max);

	void SetTurn(double angle, Vect3d axis);	// axis должен быть единичным векторм
	void AddSmallTurn(Vect3d dt);
	void AddTurn(Vect3d t);
	void Rotate(Vect3d t, Vect3d v);

	Vector3D& operator +=(Vect3d r)				{ x += r.x;	y += r.y; z += r.z; return *this; }
	Vector3D& operator -=(Vect3d r)				{ x -= r.x;	y -= r.y; z -= r.z; return *this; }
	Vector3D& operator *=(double k)				{ x *= k;	y *= k;   z *= k;	return *this; }
	Vector3D& operator /=(double k)				{ x /= k;	y /= k;   z /= k;	return *this; }

	Vector3D operator -()				const	{ Vector3D r(-x, -y, -z); return r; }
	Vector3D GetNormalized()			const	{ Vector3D v(*this); return v /= Abs(); }

	BOOL operator ==(Vect3d r)			const	{ return x == r.x && y == r.y && z == r.z; }
	BOOL operator !=(Vect3d r)			const	{ return !operator==(r); }
	BOOL IsNull()						const	{ return !(x || y || z); }

	double Sqr()						const	{ return x * x + y * y + z * z; }
	double Abs()						const;
	double GetNorm1()					const;

	void   NormalizeTurn();
	double GetTurnAngle()				const;

	Vector3D& Sum       (Vect3d r1, Vect3d r2)
	{
		x = r1.x + r2.x; y = r1.y + r2.y; z = r1.z + r2.z; return *this;
	}

	Vector3D& AddSum(Vect3d r1, Vect3d r2)
	{
		x += r1.x + r2.x; y += r1.y + r2.y; z += r1.z + r2.z; return *this;
	}

	Vector3D& SubtractSum(Vect3d r1, Vect3d r2)
	{
		x -= r1.x + r2.x; y -= r1.y + r2.y; z -= r1.z + r2.z; return *this;
	}

	Vector3D& Difference(Vect3d r1, Vect3d r2)
	{
		x = r1.x - r2.x; y = r1.y - r2.y; z = r1.z - r2.z; return *this;
	}

	Vector3D& AddDifference(Vect3d r1, Vect3d r2)
	{
		x += r1.x - r2.x; y += r1.y - r2.y; z += r1.z - r2.z; return *this;
	}

	Vector3D& Product   (Vect3d r0,  double k)
	{
		x = r0.x * k; y = r0.y * k; z = r0.z * k; return *this;
	}

	Vector3D& Product   (double k, Vect3d r0)
	{
		x = r0.x * k; y = r0.y * k; z = r0.z * k; return *this;
	}

	Vector3D& AddProduct(Vect3d r0,  double k)
	{
		x += r0.x * k; y += r0.y * k; z += r0.z * k; return *this;
	}

	Vector3D& AddProduct(double k, Vect3d r0)
	{
		x += r0.x * k; y += r0.y * k; z += r0.z * k; return *this;
	}

	Vector3D& SubtractProduct(Vect3d r0,  double k)
	{
		x -= r0.x * k; y -= r0.y * k; z -= r0.z * k; return *this;
	}

	Vector3D& SubtractProduct(double k, Vect3d r0)
	{
		x -= r0.x * k; y -= r0.y * k; z -= r0.z * k; return *this;
	}

	Vector3D& Quotient  (Vect3d r0,  double k)
	{
		x = r0.x / k; y = r0.y / k; z = r0.z / k; return *this;
	}

	Vector3D& VectorProduct   (Vect3d r1, Vect3d r2)
	{
		x = r1.y * r2.z - r1.z * r2.y; y = r1.z * r2.x - r1.x * r2.z; 
		z = r1.x * r2.y - r1.y * r2.x; return *this;
	}

	Vector3D& AddVectorProduct(Vect3d r1, Vect3d r2)
	{
		x += r1.y * r2.z - r1.z * r2.y; y += r1.z * r2.x - r1.x * r2.z;
		z += r1.x * r2.y - r1.y * r2.x;	return *this;
	}
	
	double x, y, z;
};

//_________________________________________________________________________________________

static const Vector3D VECT3D_NULL(0, 0, 0);		// нулевой вектор

//_________________________________________________________________________________________

inline Vector3D operator +(Vect3d r1, Vect3d r2)	{ Vector3D rr(r1);	return rr += r2;	}
inline Vector3D operator -(Vect3d r1, Vect3d r2)	{ Vector3D rr(r1);	return rr -= r2;	}
inline Vector3D operator *(Vect3d r,  double k)		{ Vector3D rr(r);	return rr *= k;		}
inline Vector3D operator *(double k,  Vect3d r)		{					return r * k;		}
inline Vector3D operator /(Vect3d r,  double k)		{ Vector3D rr(r);	return rr /= k;		}

inline Vector3D operator %(Vect3d r1, Vect3d r2)	// векторное произведение
{
	return Vector3D( r1.y * r2.z - r1.z * r2.y, r1.z * r2.x - r1.x * r2.z, r1.x * r2.y - r1.y * r2.x);
}

//_________________________________________________________________________________________
//
//	скал€рное произведение вектора и тензора
//	с диагональными компонентами

inline Vector3D operator &(Vect3d r1, Vect3d r2)
{ 
    return Vector3D(r1.x * r2.x, r1.y * r2.y, r1.z * r2.z);
}

//_________________________________________________________________________________________

inline double operator *(Vect3d r1, Vect3d r2)		// скал€рное произведение
{
	return r1.x * r2.x + r1.y * r2.y + r1.z * r2.z;
}

//_________________________________________________________________________________________

#ifdef _CRT_SECURE_NO_WARNINGS
_CRT_SECURE_NO_WARNINGS
#endif

#endif // ___Vector3D_h___
