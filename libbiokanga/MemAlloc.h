#pragma once
#include "./commdefs.h"

class CMemAlloc
{
	unsigned long m_LastErr;
	int m_NumEls;
	int m_ElSize;
	HANDLE m_hHeap;
	void *m_pAllocd;
public:
	CMemAlloc(void);
	CMemAlloc(int ElSize,int NumEls,DWORD HeapFlags=0);
	~CMemAlloc(void);
	void *Create(int ElSize,int NumEls,DWORD HeapFlags=0);
	void *Realloc(int NumEls,DWORD HeapFlags=0);
	void *pAlloc(void) { return(m_pAllocd); };
	DWORD LastError(void) { return(m_LastErr); };
};
