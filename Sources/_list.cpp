//  _List.cpp
/******************************************************************************

      ����� _List<Type>

        ������ �.

           2001

      ���������� 2017

******************************************************************************/

#include "_list.h"

class _ListNode
{
private:
    _ListNode* next;
    void* data;

public:
    _ListNode(void* const d, _ListNode* n = 0)
        :next(n), data(d)                       {}
 //   ~_ListNode()                                { delete next; } 

    void*       Data()                          { return data; }
    _ListNode*& Next()                          { return next; }
};

//_____________________________________________________________________________

//_ListHeader::~_ListHeader() { delete head; } (������� ����� ����������� �������)

_ListHeader::~_ListHeader()
{
    int more = Restore();
    while (more) { Delete(more); }
}

BOOL _ListHeader::Restore()
{ 
    cursor = &head;

    return head != 0; 
}

void _ListHeader::Insert(void* const d)
{
    head = new _ListNode(d, head);
    
    count++;
}

void* _ListHeader::Iterate(BOOL& More)
{ 
    if (*cursor == 0) { return 0; }

    _ListNode* _cursor_prev = *cursor;

    cursor = &(*cursor)->Next();
    More = *cursor != 0;

    return _cursor_prev->Data();
}

BOOL _ListHeader::Move(_ListHeader* const to) 
{
    if (*cursor == 0) { return 0; }

    _ListNode* n = *cursor;       // ��������� �� ������� �����

    *cursor = (*cursor)->Next();  // ���������� ����������� � ���������� �� ������� �����
    count--;                      // (������� ����� ����������)

    n->Next() = to->head;
    to->head = n;
    to->count++;

    return *cursor != 0;
}

void* _ListHeader::Delete(BOOL& More) 
{ 
    if (*cursor == 0) { return 0; }

    _ListNode* n = *cursor;       // ��������� �� ������� �����
    void* d = n->Data();

    *cursor = (*cursor)->Next();  // ���������� ����������� � ���������� �� ������� �����
    count--;                      // (������� ����� ����������)
    
    More = *cursor != 0;
    n->Next() = 0;

    delete n;
    return d;
}

void* _ListHeader::GetData() const
{
    return (*cursor)->Data();
}

//_____________________________________________________________________________
