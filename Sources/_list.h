//  _List.h
/******************************************************************************

      ����� _List<Type> ������ ����������

         ������ �.

           2001

        ���������� 2017

******************************************************************************/

#ifndef _ListTemplate_h_
#define _ListTemplate_h_

//#include "def.h"
//#include <iostream>
//using namespace std;

#define BOOL  int

class _ListNode;

class _ListHeader
{
protected:
    _ListHeader() : head(0), cursor(&head), count(0)    {}
    ~_ListHeader();

    void  Insert(void* const d);                //  �������� ������� � ������ ������
    void* Iterate(BOOL& More);                  //  ��������� ������� ��������� ������� ������; ���������� ��������� �� ������ ������� �������
    BOOL  Move(_ListHeader* const to);          //  ����������� ������� �� ������ ������ � ������
    BOOL  Restore();                            //  ��������� ������� ������ ������� ������; ���������� (1), ���� ������ �� ������
    void* Delete(BOOL& More);                   //  ������� ������� ������� �� ������

    void* GetData()     const;                  //  ���������� ��������� �� ������� �������
    int   GetCount()    const                           { return count; }

private:
    _ListNode* head;
    _ListNode** cursor;
    int count;
}; 

//_____________________________________________________________________________

template<class Type> class _List: private _ListHeader
{
public:
    _List(): _ListHeader()              {}
    ~_List()                            {}

    void  Insert(Type* const d)         { _ListHeader::Insert((void*)d); }
    Type* Iterate(BOOL& More)           { return (Type*)_ListHeader::Iterate(More); }
    BOOL  Move(_List<Type>* const to)   { return _ListHeader::Move(to); }
    BOOL  Restore()                     { return _ListHeader::Restore(); }
    Type* Delete(BOOL& More)            { return (Type*)_ListHeader::Delete(More); }
    void  Destroy()                     { BOOL more = Restore(); while(more) { delete Iterate(more); } }

    Type* operator->() const            { return (Type*)_ListHeader::GetData(); }
    int   GetCount() const              { return _ListHeader::GetCount(); }
};

//  �������� ������ �� ������ �� ����� �������� �������� ���� <Type>,
//  ��� �� �������� ������� ��������������� �������� Destroy()

//_____________________________________________________________________________

//_ListTemplate_h_
#endif
