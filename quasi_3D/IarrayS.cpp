//-----------------------------------------------------------------------------------
//	IarrayS.cpp
//-----------------------------------------------------------------------------------
//
//	Определения функций класса IArrayCase
//
//  В. Цаплин
//
//  21.04.2010
//
//	Изменен 18.11.2011
//	Изменен 04.06.2012
//-----------------------------------------------------------------------------------

#include "IarrayS.h"

//-----------------------------------------------------------------------------------

//int IArrayCase::IsInverseByteOrder = 1;

//-----------------------------------------------------------------------------------

void** IArrayCase::GetData(unsigned __int8* ind)
{
	void** k = data;

	unsigned __int8 ind_j = ind[3];

	if (k[ind_j] == NULL)
	{
		k[ind_j] = static_cast<void*> (new void*[256]);
		
		memset(k[ind_j], 0, sizeof(void*) * 256);
	}

	k = static_cast<void**> (k[ind_j]);

	ind_j = ind[2];

	if (k[ind_j] == NULL)
	{
		k[ind_j] = static_cast<void*> (new void*[256]);
		
		memset(k[ind_j], 0, sizeof(void*) * 256);

	}

	k = static_cast<void**> (k[ind_j]);

	return &k[ind[1]];
}

//-----------------------------------------------------------------------------------

int IArrayCase::GetCountS() const
{
	int CountS = 0;

	for (int j0 = 0; j0 < 0x100; j0++)
	{
		void** d0 = (void**)data[j0];

		if (d0 == NULL)
		{
			continue;
		}

		for (int j1 = 0; j1 < 0x100; j1++)
		{
			void** d1 = (void**)d0[j1];

			if (d1 == NULL)
			{
				continue;
			}

			for (int j2 = 0; j2 < 0x100; j2++)
			{
				void* d2 = d1[j2];

				if (d2 == NULL)
				{
					continue;
				}

				CountS += SumCount(d2);
			}
		}
	}

	return CountS;
}

//-----------------------------------------------------------------------------------

void IArrayCase::RemoveAll16()
{
	for (int j0 = 0; j0 < 0x100; j0++)
	{
		void* d0 = data[j0];

		if (d0 == NULL)
		{
			continue;
		}
		
		Remove(d0);
			
		data[j0] = NULL;
	}
}

//-----------------------------------------------------------------------------------

void IArrayCase::RemoveAll()
{
	for (int j0 = 0; j0 < 0x100; j0++)
	{
		void** d0 = (void**)data[j0];

		if (d0 == NULL)
		{
			continue;
		}

		for (int j1 = 0; j1 < 0x100; j1++)
		{
			void** d1 = (void**)d0[j1];

			if (d1 == NULL)
			{
				continue;
			}

			for (int j2 = 0; j2 < 0x100; j2++)
			{
				void* d2 = d1[j2];

				if (d2 == NULL)
				{
					continue;
				}

				Remove(d2);
			}
			
			delete [] d1;
		}

		delete [] d0;
			
		data[j0] = NULL;
	}
}

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------

