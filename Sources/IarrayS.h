//-----------------------------------------------------------------------------------
//	IarrayS.h
//-----------------------------------------------------------------------------------
//
//	Класс Array<TYPE>
//
//	Класс ArrayS<TYPE>
//
//	Класс IArray16<TYPE>
//
//	Класс IArray<TYPE>
//
//	Класс IArrayS<TYPE>
//
//  В. Цаплин
//
//  25.03.2008
//
//  Изменен 21.04.2010
//  Изменен 18.11.2011
//  Изменен 04.06.2012
//  Изменен 18.07.2013  template<class TYPE> class IArray : public IArrayCase
//                                                           ______(был private)
//  Изменен 28.03.2014  template<class TYPE> class ArrayS (был : Array)
//  Изменен 10.12.2017
//-----------------------------------------------------------------------------------

#ifndef ___ArrayS_H___
#define ___ArrayS_H___

//#include <tchar.h>

#ifndef LONG
#define LONG long
#endif

#ifndef INT
#define INT int
#endif

#ifndef NULL
#define NULL  0
#endif

template<class TYPE> class Array
{
public:
	Array()	: data(NULL),	count(0)	{}
	Array(int n) : data(NULL), count(0)	{ Create(n); }
	~Array()							{ RemoveAll(); }

	int  GetCount() const				{ return count; }
	void Create(int n);
	void RemoveAll();

	// Прямой доступ к элементам данных (может вернуть NULL)
	const TYPE* GetData() const			{ return static_cast<const TYPE*> (data); }
	TYPE* GetData()						{ return data; }

	TYPE  operator[](INT i) const		{ return data[i]; }
	TYPE& operator[](INT i)				{ return data[i]; }

	const Array& operator+=(const Array& Add)
	{
		int cm = count > Add.count ? Add.count : count;
		for (int i = 0; i < cm; i++) data[i] += Add.data[i]; return *this;
	}

	const Array& operator=(const TYPE& Equ)
	{ for (int i = 0; i < count; i++) data[i] = Equ; return *this; }

	const Array& operator+=(const TYPE& Add)
	{ for (int i = 0; i < count; i++) data[i] += Add; return *this; }

	const Array& operator*=(const TYPE& Mult)
	{ for (int i = 0; i < count; i++) data[i] *= Mult; return *this; }

	void AddProduct(const Array& Add, const TYPE Mult)
	{
		int cm = count > Add.count ? Add.count : count;
		for (int i = 0; i < cm; i++) data[i] += Add[i] * Mult;
	}

	double operator*(const Array& Mult);

protected:

	TYPE* data;		// фактический массив данных
	int count;		// число элементов
};

//---------------------------------------------------------------------------

template<class TYPE> void Array<TYPE>::RemoveAll()
{
	if (data != NULL) delete [] data;
	data = NULL;
	count = 0;
}

template<class TYPE> void Array<TYPE>::Create(int n)
{
	RemoveAll();
	if (n <= 0) return;
	data = new TYPE[n];
	count = n;
}

template<class TYPE> double Array<TYPE>:: operator*(const Array& Mult)
{
	double Prod = 0;
	int cm = count > Mult.count ? Mult.count : count;
	for (int i = 0; i < cm; i++) Prod += data[i] * Mult.data[i];
	return Prod;
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//   Простой массив
//   Двумерный массив
//   Трехмерный массив
//   Четырехмерный массив: см. пример

template<class TYPE> class ArrayS
{
public:
	ArrayS() :							Array0(),	dataS(NULL), countS(0)	{}
	ArrayS(int n1) :					Array0(),	dataS(NULL), countS(0)	{ Create(n1); }
	ArrayS(int n1, int n2) :			Array0(),	dataS(NULL), countS(0)	{ Create(n1, n2); }
	ArrayS(int n1, int n2, int n3) :	Array0(),	dataS(NULL), countS(0)	{ Create(n1, n2, n3); }
	~ArrayS()																	{ RemoveAll(); }

	//   Число элементов предыдущего ранка (ArrayS<TYPE> или TYPE)
	int  GetCount() const				{ return Array0.GetCount() ? Array0.GetCount() : countS; }

	//   Число элементов (TYPE)
	int  GetCountS() const	{ return Array0.count ? Array0.count : dataS ? countS * dataS->GetCountS() : 0; }

	//   Размерность (количество измерений) массива
	int  GetRank() const				{ return Array0.GetCount() ? 1 : dataS ? 1 + dataS->GetRank() : 0; }

	void Create(int n1)					{ RemoveAll(); Array0.Create(n1); }
	void Create(int n1, int n2)			{ CreateS(n1); for (int i = 0; i < n1; i++) dataS[i].Create(n2); }
	void Create(int n1, int n2, int n3)	{ CreateS(n1); for (int i = 0; i < n1; i++) dataS[i].Create(n2, n3); }
	void RemoveAll();

	// Прямой доступ к элементам данных (может вернуть NULL)
	const ArrayS* GetDataS() const		{ return static_cast<const ArrayS*> (dataS); }
	ArrayS* GetDataS()					{ return static_cast<ArrayS*> (dataS); }

	ArrayS& operator[](LONG i)			{ return dataS[i]; }	// для 2D & 3D массивов
																// (2D : Элемент == ArrayName[LONG][INT])
																// (3D : Элемент == ArrayName[LONG][LONG][INT])

	TYPE&   operator[](INT i)			{ return Array0.GetData()[i]; }	// Для простых массивов
																		// (Элемент == ArrayName[INT])

	void AddProduct(const Array<TYPE>& V1, const Array<TYPE>& V2);
	void Product(const ArrayS<TYPE>& M, const Array<TYPE>& V);

	const ArrayS& operator=(const TYPE& Equ);

	const ArrayS& operator*=(const TYPE& Mult);
	const ArrayS& operator+=(const TYPE& Add);
	const ArrayS& operator+=(const ArrayS<TYPE>& Add);

protected:
	void CreateS(int n1);

	Array<TYPE> Array0;
	ArrayS* dataS;	// фактический массив данных второго и более уровня
	int countS;		// число элементов второго и более уровня
};

//	Создание черырехмерного массива:
//	{
//		...
//		ArrayS<double> FourDimensionalArray(Size1, 0);
//		for (i1 = 0; i1 < Size1; i1++)	FourDimensionalArray[i1].Create(Size2, Size3, Size4);
//		...
//	}

//---------------------------------------------------------------------------

template<class TYPE> void ArrayS<TYPE>::RemoveAll()
{
	if (dataS != NULL) delete [] dataS;
	dataS = NULL;
	countS = 0;

	Array0.RemoveAll();
}

template<class TYPE> void ArrayS<TYPE>::CreateS(int n1)
{
	RemoveAll();
	if (n1 <= 0) return;
	dataS = new ArrayS[n1];
	countS = n1;
}

template<class TYPE> void ArrayS<TYPE>::AddProduct(const Array<TYPE>& V1, const Array<TYPE>& V2)
{
	int cm1 = V1.GetCount(),  cm2 = V2.GetCount();
	if (cm1 != countS)  return;
	for (int i = 0; i < cm1; i++)
	{	if (cm2 != dataS[i].GetCount())  return;
		for (int j = 0; j < cm2; j++)  dataS[i][j] += V1[i] * V2[j]; }
}

template<class TYPE> void ArrayS<TYPE>::Product(const ArrayS<TYPE>& M, const Array<TYPE>& V)
{
	if (M.Array0.GetCount() != 0 || M.dataS->Array0.GetCount() != V.GetCount())  return;
	if (Array0.GetCount() == 0 || Array0.GetCount() != M.countS)  Create(M.countS);
	for (int i = 0; i < Array0.GetCount(); i++) Array0.GetData()[i] = M.dataS[i].Array0 * V;
}

template<class TYPE> const ArrayS<TYPE>& ArrayS<TYPE>::operator=(const TYPE& Equ)
{
	if (Array0.GetCount())
	{ for (int i = 0; i < Array0.GetCount(); i++) Array0.GetData()[i] = Equ; }
	else
	{ for (long i = 0; i < countS; i++) dataS[i] = Equ; }
	return *this;
}

template<class TYPE> const ArrayS<TYPE>& ArrayS<TYPE>::operator*=(const TYPE& Mult)
{
	if (Array0.GetCount())
	{ for (int i = 0; i < Array0.GetCount(); i++) Array0.GetData()[i] *= Mult; }
	else
	{ for (long i = 0; i < countS; i++) dataS[i] *= Mult; }
	return *this;
}

template<class TYPE> const ArrayS<TYPE>& ArrayS<TYPE>::operator+=(const TYPE& Add)
{
	if (Array0.GetCount())
	{ for (int i = 0; i < Array0.GetCount(); i++) Array0.data[i] += Add; }
	else
	{ for (long i = 0; i < countS; i++) dataS[i] += Add; }
	return *this;
}

template<class TYPE> const ArrayS<TYPE>& ArrayS<TYPE>::operator+=(const ArrayS<TYPE>& Add)
{
	if (Array0.GetCount())
	{ for (int i = 0; i < Array0.GetCount(); i++) Array0.GetData()[i] += Add.Array0.GetData()[i]; }
	else if (countS == Add.countS)
	{ for (long i = 0; i < countS; i++) dataS[i] += Add.dataS[i]; }
	return *this;
}

//---------------------------------------------------------------------------

//   Бесконечный массив ([int16], [int32] элементов)

//---------------------------------------------------------------------------

#include <memory.h>

class IArrayCase
{
protected:

	IArrayCase()						{ memset(data, 0, sizeof(void*) * 256); }

	void** GetData  (unsigned __int8 * index);

	void* data[256];					// фактический массив данных

	int  GetCountS() const;
	void RemoveAll();
	void RemoveAll16();

	//static int IsInverseByteOrder;

private:

	virtual int  SumCount(void* el) const = 0;
	virtual void Remove(void* el) = 0;
};

//---------------------------------------------------------------------------

//   Бесконечный массив 16

//---------------------------------------------------------------------------

template<class TYPE> class IArray16 : private IArrayCase
{
public:
	IArray16()	: count(0)			{}
	~IArray16()						{ if (count) RemoveAll(); }

	int  GetCount() const			{ return count; }
	void RemoveAll()				{ IArrayCase::RemoveAll16(); count = 0; }

	TYPE& operator[](INT i);

protected:

	int count;						// число элементов
	int  SumCount(void* el) const	{ return 0; }
	void Remove(void* el)			{ delete [] reinterpret_cast<TYPE*> (el); }
};

//---------------------------------------------------------------------------

template<class TYPE> TYPE& IArray16<TYPE>::operator[](INT i)
{
	unsigned __int8& index = *reinterpret_cast<unsigned __int8*> (&i);

	void** ppData = &data[(&index)[1]];

	if (*ppData == NULL)
	{
		*ppData = static_cast<void*> (new TYPE[256]);
		count += 256;
	}

	return (static_cast<TYPE*> (*ppData))[index];
}

//---------------------------------------------------------------------------

//   Бесконечный массив 32

//---------------------------------------------------------------------------

template<class TYPE> class IArray : public IArrayCase
{
public:
	IArray()	: count(0)			{}
	~IArray()						{ if (count) RemoveAll(); }

	int  GetCount() const			{ return count; }
	void RemoveAll()				{ IArrayCase::RemoveAll(); count = 0; }

	TYPE& operator[](INT i);

protected:

	int count;						// число элементов
	int  SumCount(void* el) const	{ return 0; }
	void Remove(void* el)			{ delete [] reinterpret_cast<TYPE*> (el); }
};

//---------------------------------------------------------------------------

template<class TYPE> TYPE& IArray<TYPE>::operator[](INT i)
{
	unsigned __int8& index = *reinterpret_cast<unsigned __int8*> (&i);

	void** ppData = GetData(&index);

	if (*ppData == NULL)
	{
		*ppData = static_cast<void*> (new TYPE[256]);
		count += 256;
	}

	return (static_cast<TYPE*> (*ppData))[index];
}

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

//   Многомерный бесконечный массив

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------

template<class TYPE> class IArrayS : private IArray<TYPE>
{
public:
	IArrayS() :	IArray<TYPE> (), countS(0)	{}
	~IArrayS()								{ RemoveAll(); }

	//   Число элементов предыдущего ранка (ArrayS<TYPE> или TYPE)
	//int  GetCount() const					{ return count ? count	: countS; }
	int  GetCount() const					{ return IArray<TYPE>::GetCount() ? IArray<TYPE>::GetCount() : countS; }

	//   Число элементов (TYPE)
	int  GetCountS() const;

	void RemoveAll()						{ IArrayCase::RemoveAll(); countS = 0; }
	//void RemoveAll()						{ IArray<TYPE>::RemoveAll(); countS = 0; }

	IArrayS& operator[](LONG i);
	TYPE&    operator[](INT i);

protected:

	int countS;		// число элементов

	//static IArrayS dummy;
	static TYPE dummy;

	int  SumCount(void* el) const;
	void Remove(void* el);
};

//---------------------------------------------------------------------------

template<class TYPE> TYPE IArrayS<TYPE>::dummy;

//---------------------------------------------------------------------------

template<class TYPE> TYPE& IArrayS<TYPE>::operator[](INT i)
{
	if (countS)
	{
		return dummy;
	}

	return IArray<TYPE>::operator[](i);
}

//---------------------------------------------------------------------------

template<class TYPE> IArrayS<TYPE>& IArrayS<TYPE>::operator[](LONG i)
{
	if (IArray<TYPE>::GetCount())
	{
		return *this;
	}

	unsigned __int8& index = *reinterpret_cast<unsigned __int8*> (&i);

	void** ppData = IArrayCase::GetData(&index);

	if (*ppData == NULL)
	{
		*ppData = static_cast<void*> (new IArrayS[256]);
		countS += 256;
	}

	return (static_cast<IArrayS*>(*ppData))[/*IsInverseByteOrder == 1 ? */index/* : *(&index + 3)*/];
}

//---------------------------------------------------------------------------

template<class TYPE> int IArrayS<TYPE>::GetCountS() const
{
	if (IArray<TYPE>::GetCount())
	{
		return IArray<TYPE>::GetCount();
	}

	if (!countS)
	{
		return 0;
	}

	return IArrayCase::GetCountS();
}

//---------------------------------------------------------------------------

template<class TYPE> int IArrayS<TYPE>::SumCount(void* el) const
{
	int CountS = 0;

	IArrayS* El = static_cast<IArrayS*> (el);

	for (int j3 = 0; j3 < 0x100; j3++)
	{
		CountS += El[j3].GetCountS();
	}

	return CountS;
}

//---------------------------------------------------------------------------

template<class TYPE> void IArrayS<TYPE>::Remove(void* el)
{
	if (IArray<TYPE>::GetCount())
	{
		delete [] reinterpret_cast<TYPE*> (el);
	}
	else if (countS)
	{
		delete [] reinterpret_cast<IArrayS*> (el);
	}
}

//---------------------------------------------------------------------------

//template <class type> inline void swap(type &a, type &b) { type c = a; a = b; b = c; }

//---------------------------------------------------------------------------

#endif //___ArrayS_H___

//---------------------------------------------------------------------------
