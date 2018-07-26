//  _List.h
/******************************************************************************

      Класс _List<Type> списка переменных

         Цаплин В.

           2001

        Исправлено 2017

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

    void  Insert(void* const d);                //  вставить элемент в начало списка
    void* Iterate(BOOL& More);                  //  назначить текущим следующий элемент списка; возвращает указатель на бывший текущий элемент
    BOOL  Move(_ListHeader* const to);          //  переместить элемент из одного списка в другой
    BOOL  Restore();                            //  назначить текущим первый элемент списка; возвращает (1), если список не пустой
    void* Delete(BOOL& More);                   //  удалить текущий элемент из списка

    void* GetData()     const;                  //  возвращает указатель на текущий элемент
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

//  удаление списка не влечет за собой удаление объектов типа <Type>,
//  для их удаления следует воспользоваться функцией Destroy()

//_____________________________________________________________________________

//_ListTemplate_h_
#endif
