//	util.cpp

#include "util.h"
#include <stdio.h>

//------------------------------------------------------------------------..

int is_inverse_byte_order()
{
	short x = 1 ;
	char *p = reinterpret_cast<char*> (&x);
	return 0 != (int)(*p);
}

//------------------------------------------------------------------------..

int InputInt(FILE* file, int& err)
{
	char inp_s[300];
	int inp_i;

	if (fscanf(file, "\n%[^\n]\n", inp_s) != 1)
	{
		err = 1;
	}

	if (sscanf(inp_s, "%i", &inp_i) != 1)
	{
		err = 1;
	}

	return inp_i;
}

//_________________________________________________________________________________

double InputDouble(FILE* file, int& err)
{
	char inp_s[300];

	if (fscanf(file, "\n%[^\n]\n", inp_s) != 1)
	{
		err = 1;
	}

	return atof(inp_s);
}

//_________________________________________________________________________________

char* dos_rus(const char* txt, char* buf)
{
    unsigned char* p = (unsigned char*)txt;
    unsigned char* pb = (unsigned char*)buf;

    for (; p < (unsigned char*)(txt + strlen(txt)); p++, pb++)
    {
        if (*p >=192 && *p < 240)       { *pb = *p - 64; }
        else if (*p >= 240)             { *pb = *p - 16; }
        else if (*p == 168)             { *pb = 240; }
        else if (*p == 184)             { *pb = 241; }
        else                            { *pb = *p; }
    }

    *pb = 0;
    return buf;
}

//_________________________________________________________________________________

/*void print(const char* format)
{
    char str[256];
    printf(dos_rus(format, str));
}

void print(const char* format, double d)
{
    char str[256];
    printf(dos_rus(format, str), d);
}

void print(const char* format, int i)
{
    char str[256];
    printf(dos_rus(format, str), i);
}

void print(const char* format, ...)
{
    char str[256];
    const char* a=strchr(format, '%');
    if (a==0) { printf(dos_rus(format, str)); }
    else if (*(a+1)=='f' || *(a+1)=='e' || *(a+1)=='.')
    {
        double d=*(double*)(&format+1);
        printf(dos_rus(format, str), d);
    }
    else if (*(a+1)=='i' || *(a+1)=='0')
    {
        int i=*(int*)(&format+1);
        printf(dos_rus(format, str), i);
    }
    else if (*(a+1)=='u')
    {
        unsigned int u=*(unsigned int*)(&format+1);
        printf(dos_rus(format, str), u);
    }
    else  { printf(dos_rus(format, str)); }
}*/

void print(const char* format, ...)
{
    //char str[256];  dos_rus(format, str);
    //char* a=str;
    //char* b=str;
    //// ищем ближайший '%', за которым сразу не следует еще один '%':
    //do { b=strchr(b, '%'); } while (b && *((++b)++)=='%');
    //if (b==0) { printf(str); return; }
    //b-=2;
    //char *c, *argum = (char*)(&format+1);

    //do  {
    //    *b='%'; c=b+1; // восстанавливаем временно обнуленный символ (если это не первый цикл)
    //    // ищем следующий '%', за которым сразу не следует еще один '%':
    //    do { c=strchr(c, '%'); } while (c && *((++c)++)=='%');
    //    if (c) { c-=2; *c=0; } // временно укорачиваем строку format

    //    if (*(b+1)=='f' || *(b+1)=='e' || *(b+1)=='.')
    //    {
    //        double d=*(double*)argum; argum+=8;//sizeof(double)
    //        printf(a, d);
    //    }
    //    else if (*(b+1)=='i' || *(b+1)=='0' || *(b+1)=='u')
    //    {
    //        int i=*(int*)argum; argum+=4;//sizeof(int)
    //        printf(a, i);
    //    }
    //    else if (*(b+1)=='L')
    //    {
    //        __int64 L=*(__int64*)argum; argum+=8;//sizeof(__int64)
    //        char s1[50];
    //        _i64toa(L, s1, 10);
    //        *(b+1) = 's';
    //        printf(a, s1);
    //    }
    //    else if (*(b+1)=='s')
    //    {
    //        char *s=*(char**)argum; argum+=4;//sizeof(char*)
    //        char s2[256];
    //        printf(a, dos_rus(s, s2));
    //    }
    //    else
    //    { printf("Error \" print\"\n"); return; }

    //    a=b=c; // оставшаяся невыведенная часть строки format
    //} while (a);
}
