//-----------------------------------------------------------------------------------
//
//	Mathematical utilities
//
//	Anton M. Krivtsov
//
//	11.04.2001
//
//-----------------------------------------------------------------------------------

#ifndef ___UTIL_H___
#define ___UTIL_H___
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//#define __USE_STD_IOSTREAM

//#include <tchar.h>
//#include <sstream.h>

//using namespace std;
//-----------------------------------------------------------------------------------

//static const double DEG = PI / 180;

//-----------------------------------------------------------------------------------

int is_inverse_byte_order();

//-----------------------------------------------------------------------------------
//	RAND_MAX = 32767 = 2^15 - 1
//-----------------------------------------------------------------------------------

inline double rand(const double max)			//	Random double value from 0 to max
{
	return max * (double)rand() / RAND_MAX;
}

//-----------------------------------------------------------------------------------

//inline double rand(const double min, const double max)	//	Random double value from min to max
//{
//	return min + (max - min) * (double)rand() / RAND_MAX;
//}

//-----------------------------------------------------------------------------------

inline int rand(int max)		//	Random integer value from 0 to (max - 1)
{
	return max * rand() >> 15;	//	x >> 15 = x / ((int)RAND_MAX + 1)
}

//-----------------------------------------------------------------------------------

template <class T> inline void swap(T& x1, T& x2) { T x = x1; x1 = x2; x2 = x; }

//-----------------------------------------------------------------------------------

#ifndef sqr
#define sqr square
template <class T> inline T sqr(const T& x) { return x * x; }
#endif

//-----------------------------------------------------------------------------------

/*template <class T> string str_form(const T& x)
{
	stringstream buf;
	buf << x;
	return buf.str();
}*/

//-----------------------------------------------------------------------------------

template <class type> type max_num(const type a, const type b) { return a>b? a:b; }
template <class type> type min_num(const type a, const type b) { return a<b? a:b; }

template <class type> type max_num(const type a, const type b, const type c) { return max_num(a, max_num(b, c) ); }
template <class type> type min_num(const type a, const type b, const type c) { return min_num(a, min_num(b, c) ); }

template <class type> type min_max_num(const type x, const type a, const type b) { return max_num( min_num(x, b), a); }

//-----------------------------------------------------------------------------------

inline void swap_float(float &x)
{
  char *p = reinterpret_cast<char*>(&x);
  char q[4] = {p[3], p[2], p[1], p[0]};
  memcpy((void*)(&x), (void*)q, 4);
}

//-----------------------------------------------------------------------------------

inline void swap_double(double &x)
{
  char *p = reinterpret_cast<char*>(&x);
  char q[8] = {p[7], p[6], p[5], p[4], p[3], p[2], p[1], p[0]};
  memcpy((void*)(&x), (void*)q, 8);
}

//-----------------------------------------------------------------------------------

int    InputInt   (FILE* file, int& err);
double InputDouble(FILE* file, int& err);

//-----------------------------------------------------------------------------------

char* dos_rus(const char* txt, char* buf);
//void print(const char* format);
//void print(const char* format, double d);
//void print(const char* format, int i);
void print(const char* format, ...);

//-----------------------------------------------------------------------------------
#endif //___UTIL_H___

