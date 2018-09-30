
//	Тензорные линейные операции

//	tlinear3D.h

//	19.02.2008

//  Изменен 05.04.2010
//  Изменен 14.09.2011
//  Изменен 21.12.2011

//___________________________________________________________________________________________________

#include "vector3d.h"
//#include <tchar.h>
#include <memory.h>

//___________________________________________________________________________________________________

#ifndef ___TLINEAR
#define ___TLINEAR

//___________________________________________________________________________________________________
//
//
//	Класс трехмерных поворотов
//

class Turn
{
public:

	Turn(): w(1), x(0), y(0), z(0)				{}
	Turn(const Turn& t)							{ w = t.w; x = t.x; y = t.y; z = t.z; }
	Turn(double w, double x = 0, double y = 0, double z = 0)
												{ Set(w, x, y, z); }
	Turn(Vect3d v)								{ SetTurnVector(v); }

	void SetUnit()								{ w = 1; x = 0; y = 0; z = 0; }
	void Set(Vect3d r)							{ x = r.x; y = r.y; z = r.z; }
	void Set(double w, double x, double y, double z)
	{
		this->w = w; this->x = x; this->y = y; this->z = z;
	}

	void SetTurnVector(Vect3d v);
	void SetRand();

	// операции с кватернионами:

	Turn& operator +=(const Turn& t)			{ w += t.w; x += t.x; y += t.y; z += t.z; return *this; }
	Turn& operator -=(const Turn& t)			{ w -= t.w; x -= t.x; y -= t.y; z -= t.z; return *this; }
	Turn& operator *=(double k)					{ w *= k;   x *= k;   y *= k;   z *= k;   return *this; }
	Turn& operator /=(double k)					{ w /= k;   x /= k;   y /= k;   z /= k;   return *this; }

	Turn operator -()					const	{ return Turn(-w, -x, -y, -z); }
	Turn operator !()					const	{ return Turn(w, -x, -y, -z); }
	Turn operator ~()					const	{ return Turn(w, -x, -y, -z) /= NormSqr(); }

	BOOL operator ==(const Turn& t)		const	{ return w == t.w && x == t.x && y == t.y && z == t.z; }
	BOOL operator !=(const Turn& t)		const	{ return !operator==(t); }

	Turn& SetProduct(const Turn& t1, const Turn& t2);
	Turn& SetProduct(const Turn& t1, Vect3d      v2);
	Turn& SetProduct(Vect3d      v1, const Turn& t2);

	double NormSqr()					const	{ return w * w + x * x + y * y + z * z; }
	double Norm()						const;

	// операции с поворотами

	Turn&    Determine(double angle, Vect3d axis);		// axis должен быть единичным вектром!

	Vector3D	ConvertToVector3D()		const	{ return Vector3D(x, y, z); }

	Vector3D	GetTurnVector()			const;	// возвращает нормализованный вектор поворота

	double		GetAngle()				const;	// возвращает от 0 до 2 * PI (для Normalized1 - от 0 до PI)
	Turn		GetNormalized1()		const;
	Turn		GetNormalized2()		const;

	Turn& SetComposed(const Turn& t1, const Turn& t2);
	Turn& SetComposed(const Turn& t1, Vect3d      v2);	// вектор малого поворота v2 = angle * axis
	Turn& SetComposed(Vect3d      v1, const Turn& t2);	// вектор малого поворота v1 = angle * axis

	Turn& Normalize1();
	Turn& Normalize2(double* pScal = 0);

	double w, x, y, z;
private:
};

//___________________________________________________________________________________________________
//
//	Неопределенный поворот

static const Turn TURN_IND = Turn(0,0,0,0);

//___________________________________________________________________________________________________

inline Turn operator +(const Turn& t1, const Turn& t2)	{ Turn tt(t1); return tt += t2;	}
inline Turn operator -(const Turn& t1, const Turn& t2)	{ Turn tt(t1); return tt -= t2;	}
inline Turn operator *(const Turn& t,  double       k)	{ Turn tt(t);  return tt *= k;	}
inline Turn operator *(double      k,  const  Turn& t)	{              return t * k;	}
inline Turn operator /(const Turn& t,  double       k)	{ Turn tt(t);  return tt /= k;	}

inline Turn operator %(const Turn& t1, Vect3d     v2)	{ return Turn().SetProduct(t1, v2); }
inline Turn operator %(Vect3d     v1, const Turn& t2)	{ return Turn().SetProduct(v1, t2); }
inline Turn operator *(const Turn& t1, const Turn& t2)	{ return Turn().SetProduct(t1, t2); }

//___________________________________________________________________________________________________
//
//	Операции с векторным базисом
//

class Base3DLinearObject
{
public:
	Base3DLinearObject();
	Base3DLinearObject(Vect3d first, Vect3d second, Vect3d third);

	Vector3D& GetFirst()					{ return x; }
	Vector3D& GetSecond()					{ return y; }
	Vector3D& GetThird()					{ return z; }

	Base3DLinearObject GetReciproca()	const;
	Vector3D ThisToNormal(Vect3d v)		const;
	Vector3D NormalToThis(Vect3d v)		const;

	Vector3D x, y, z;
};

//___________________________________________________________________________________________________

Vector3D Convert(const Base3DLinearObject& basis_dst, const Base3DLinearObject& basis_src, Vect3d v);

//___________________________________________________________________________________________________

static const Vector3D _3X3UNITM[3] = { Vector3D(1,0,0), Vector3D(0,1,0), Vector3D(0,0,1) };	// 3x3 единичная матрица
static const Vector3D _3X3NULLM[3] = { Vector3D(0,0,0), Vector3D(0,0,0), Vector3D(0,0,0) };	// 3x3 нулевая матрица

//___________________________________________________________________________________________________
//
//	Класс трехмерного
//	тензора второго ранга
//

#define Tens const Tensor&

class Tensor
{
public:
	typedef Vector3D MATRIX_3X3[3];
	typedef const MATRIX_3X3 CMATRIX_3X3;
	typedef Vector3D *PMATRIX_3X3;

	Tensor()									{}
	Tensor(Vect3d x, Vect3d y, Vect3d z)		{ Set(x, y, z); }
	Tensor(CMATRIX_3X3 m)						{ Set(m); }
	Tensor(const Turn& t)						{ Set(t); }
	Tensor(Vect3d v)							{ Turn t(v); Set(t); }
	Tensor(double a)							{ x.Set(a,a,a); y.Set(a,a,a); z.Set(a,a,a); }
	//Tensor(const Tensor& t)						{ *this = t; }

	void Set(const Turn& t);
	void Set(CMATRIX_3X3 m)						{ Set(m[0], m[1], m[2]); }
	void Set(Vect3d x, Vect3d y, Vect3d z)		{ this->x = x; this->y = y; this->z = z; }

	void SetRand(double norm_mean);
	void SetRand(double min1, double min2, double min3, double max1, double max2, double max3);

	void SetRand1()								{ SetRand(1.0); *this /= GetNorm2(); }
	void SetRand2();
	void SetRand3();

	Turn   GetTurn()					const;

	Tensor& operator +=(Tens t)					{ x += t.x; y += t.y; z += t.z; return *this; }
	Tensor& operator -=(Tens t)					{ x -= t.x; y -= t.y; z -= t.z; return *this; }
	Tensor& operator *=(double k)				{ x *= k;   y *= k;   z *= k;   return *this; }
	Tensor& operator /=(double k)				{ x /= k;   y /= k;   z /= k;   return *this; }

	Tensor  GetReciproca(double* _Det = NULL)	const;
	Tensor& SetReciproca(Tens t, double* _Det = NULL);

	Tensor  GetTransposed()				const;
	Tensor& SetTransposed(Tens t);

	Tensor operator -()					const	{ return Tensor(-x, -y, -z); }

	Tensor GetNormalized(double* _a1 = NULL,
                         double* _a2 = NULL, double* _a3 = NULL)	const;

	Vector3D GetVector()				const	{ return Vector3D(y.z - z.y, z.x - x.z, x.y - y.x); }

	double GetDeterminant()				const	{ return x % y * z; }
	double GetDistortion1()				const
	{
		return (Tensor(_3X3UNITM) -= _Transpose(*this)).GetNorm1();
	}

	double GetDistortion2()				const
	{
		return (Tensor(_3X3UNITM) -= _Transpose(*this)).GetNorm2();
	}

	double GetNorm1()					const;
	double GetNorm2()					const;
	double GetNormSqr()					const	{ return (x.Sqr() + y.Sqr() + z.Sqr()) / 3; }
	double GetSkewness1()				const	{ return GetVector().GetNorm1(); }
	double GetSkewness2()				const	{ return GetVector().Abs(); }
	double GetTrace()					const	{ return x.x + y.y + z.z; }

	Tensor Pow(int n)					const;

//	квадратный корень симметричного тензора
	Tensor Sqrt()						const;

//	квадратный корень тензора
	Tensor SqrtPlus(int* _steps = NULL, double* _det = NULL)	const;

	Tensor& Sum       (Tens t1, Tens t2)
	{
		x.Sum(t1.x, t2.x);
		y.Sum(t1.y, t2.y);
		z.Sum(t1.z, t2.z);

		return *this;
	}

	Tensor& Difference(Tens t1, Tens t2)
	{
		x.Difference(t1.x, t2.x);
		y.Difference(t1.y, t2.y);
		z.Difference(t1.z, t2.z);

		return *this;
	}

	Tensor& Product   (Tens t0, double k)
	{
		x.Product(t0.x, k);
		y.Product(t0.y, k);
		z.Product(t0.z, k);

		return *this;
	}

	Tensor& Product   (double k, Tens t0)
	{
		x.Product(t0.x, k);
		y.Product(t0.y, k);
		z.Product(t0.z, k);

		return *this;
	}

	Tensor& Quotient  (Tens t0, double k)
	{
		x.Quotient(t0.x, k);
		y.Quotient(t0.y, k);
		z.Quotient(t0.z, k);

		return *this;
	}

	Tensor& ScalarProduct	 (Tens t1, Tens t2);
	Tensor& ScalarProductTr	 (Tens t1, Tens t2);	// == ScalarProduct(t1, t2.GetTransposed())

	Tensor& VectorProduct	 (Tens t1, Vect3d v2)
	{
		Tensor t;

		t.x.VectorProduct(t1.x, v2);
		t.y.VectorProduct(t1.y, v2);
		t.z.VectorProduct(t1.z, v2);

		memcpy(this, &t, sizeof(Tensor));

		return *this;
	}

	Tensor& VectorProduct	 (Vect3d v1, Tens t2)
	{
		Tensor t;

		t.x.Product        (v1.y, t2.z);
		t.x.SubtractProduct(v1.z, t2.y);

		t.y.Product        (v1.z, t2.x);
		t.y.SubtractProduct(v1.x, t2.z);

		t.z.Product        (v1.x, t2.y);
		t.z.SubtractProduct(v1.y, t2.x);

		memcpy(this, &t, sizeof(Tensor));

		return *this;
	}

	Tensor& TensorProduct	 (Vect3d v1, Vect3d v2)
	{
		Tensor t;

		t.x.Product(v2, v1.x);
		t.y.Product(v2, v1.y);
		t.z.Product(v2, v1.z);

		memcpy(this, &t, sizeof(Tensor));

		return *this;
	}

	Tensor& AddProduct		 (Vect3d v1, Vect3d v2, double a)
	{
		x.AddProduct(v2, v1.x * a);
		y.AddProduct(v2, v1.y * a);
		z.AddProduct(v2, v1.z * a);

		return *this;
	}

	Tensor& SubtractProduct	 (Vect3d v1, Vect3d v2, double a)
	{
		x.SubtractProduct(v2, v1.x * a);
		y.SubtractProduct(v2, v1.y * a);
		z.SubtractProduct(v2, v1.z * a);

		return *this;
	}

	Vector3D x, y, z;

private:

	Vector3D _SProd(Vect3d r)				const	{ return Vector3D(x * r, y * r, z * r); }

	Tensor   _Transpose(Tens t2)	const
	{
		return Tensor(t2._SProd(x), t2._SProd(y), t2._SProd(z));
	}
};

#define GetDistortion GetDistortion1
#define GetNorm GetNorm1
#define GetSkewness GetSkewness1

//___________________________________________________________________________________________________

static const Tensor UNIT_T2R = _3X3UNITM;	// единичный тензор второго ранга
static const Tensor NULL_T2R = _3X3NULLM;	// нулевой тензор второго ранга

//___________________________________________________________________________________________________

inline Tensor operator +(Tens t1, Tens t2)			{ Tensor tt(t1); return tt += t2; }
inline Tensor operator -(Tens t1, Tens t2)			{ Tensor tt(t1); return tt -= t2; }
inline Tensor operator *(Tens t,  double k)			{ Tensor tt(t);  return tt *= k;  }
inline Tensor operator *(double k,  Tens t)			{                return t * k;    }
inline Tensor operator /(Tens t,  double k)			{ Tensor tt(t);  return tt /= k;  }

inline Vector3D operator *(Tens t,  Vect3d r)
{
	return Vector3D(t.x * r, t.y * r, t.z * r);
}

inline Vector3D operator *(Vect3d r,  Tens t)
{
	return r.x * t.x + r.y * t.y + r.z * t.z;
}

inline Tensor operator *(Tens t1, Tens t2)
{
	return Tensor(t1.x * t2, t1.y * t2, t1.z * t2);
}

//___________________________________________________________________________________________________
//
//	tensor1 / tensor2 == tensor1._Transpose(tensor2) == tensor1 * tensor2.GetTransposed()
//
inline Tensor operator /(Tens t1, Tens t2)
{
	return Tensor(t2 * t1.x, t2 * t1.y, t2 * t1.z);
}

//___________________________________________________________________________________________________
//
//	Тензорное произведение двух тензоров
//
inline Tensor operator ,(Vect3d r1, Vect3d r2)
{
	return Tensor(r2 * r1.x, r2 * r1.y, r2 * r1.z);
}

//___________________________________________________________________________________________________
//
//	Векторное произведение тензора и вектора
//

inline Tensor operator %(Tens t,  Vect3d r)
{
	return Tensor(t.x % r, t.y % r, t.z % r);
}

inline Tensor operator %(Vect3d r,  Tens t)
{
	return Tensor(r.y * t.z - r.z * t.y,  r.z * t.x - r.x * t.z,  r.x * t.y - r.y * t.x);
}

//___________________________________________________________________________________________________

inline Tensor operator &(Tens t,  Vect3d r)
{
	return Tensor(t.x & r, t.y & r, t.z & r);
}

inline Tensor operator &(Vect3d r,  Tens t)
{
	return Tensor(t.x * r.x, t.y * r.y, t.z * r.z);
}

//_________________________________________________________________________________________
// опртимизация скорости:

inline Vector3D& ScalarProduct(Vector3D& v, Tens t1, Vect3d v2)
{
	double x, y;

	x   = t1.x * v2;
	y   = t1.y * v2;
	v.z = t1.z * v2;

	v.x = x;
	v.y = y;

	return v;
}

inline Vector3D& ScalarProduct(Vector3D& v, Vect3d v1, Tens t2)
{
	Vector3D v_t;

	v_t.Product   (v1.x, t2.x);
	v_t.AddProduct(v1.y, t2.y);
	v_t.AddProduct(v1.z, t2.z);

	v.x = v_t.x;
	v.y = v_t.y;
	v.z = v_t.z;

	return v;
}

inline Tensor& Tensor::ScalarProduct  (Tens t1, Tens t2)
{
	Tensor t;

	::ScalarProduct(t.x, t1.x, t2);
	::ScalarProduct(t.y, t1.y, t2);
	::ScalarProduct(t.z, t1.z, t2);

	memcpy(this, &t, sizeof(Tensor));

	return *this;
}

inline Tensor& Tensor::ScalarProductTr(Tens t1, Tens t2)	// == ScalarProduct(t1, t2.GetTransposed())
{
	Tensor t;

	::ScalarProduct(t.x, t2, t1.x);
	::ScalarProduct(t.y, t2, t1.y);
	::ScalarProduct(t.z, t2, t1.z);

	memcpy(this, &t, sizeof(Tensor));

	return *this;
}

//___________________________________________________________________________________________________
//
//	Главные значения тензора второго ранга;
//	Тензор главных осей тензора второго ранга

//  три действительных числа

void  MainValueT2R(double& a1, double& a2, double& a3, Tens t);
void  MainValueT2R(double& a1, double& a2, double& a3, Tensor& m, Tens t);// t - симметричный

//  одно действительное число

void  MainValueT2R(double& a1, Tens t);
void  MainValueT2R(double& a1, Vector3D& m, Tens t);

//___________________________________________________________________________________________________

inline double fabs(const Tensor& T)
{
	return T.GetNorm();
}

//___________________________________________________________________________________________________

#endif
