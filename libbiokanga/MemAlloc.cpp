/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CMemAlloc::CMemAlloc(void)
{
m_hHeap = NULL;
m_NumEls = 0;
m_ElSize = 0;
m_pAllocd = NULL;
m_LastErr = 0;
}


CMemAlloc::CMemAlloc(int ElSize,int NumEls,DWORD HeapFlags)
{
m_hHeap = NULL;
m_NumEls = 0;
m_ElSize = 0;
m_pAllocd = NULL;
m_LastErr = 0;
Create(ElSize,NumEls,HeapFlags);
}

CMemAlloc::~CMemAlloc(void)
{
if(m_hHeap != NULL)
	HeapDestroy(m_hHeap);
}

void *
CMemAlloc::Create(int ElSize,int NumEls,DWORD HeapFlags)
{
void *pAllocd;
m_LastErr = 0;
if(ElSize < 1)			// make the initial allocation worth the exercise
	ElSize = 1;
if((NumEls * ElSize) < 1024)
	NumEls = (1024/ElSize) + 1;
if(m_hHeap == NULL)		// create heap, or reuse if one already created
	{
	if((m_hHeap = HeapCreate(HeapFlags,ElSize * NumEls,0))==NULL)
		{
		m_LastErr = GetLastError();
		return(NULL);
		}
	pAllocd = HeapAlloc(m_hHeap,HeapFlags,ElSize*NumEls);
	}
else
	pAllocd = HeapReAlloc(m_hHeap,0,m_pAllocd,ElSize*NumEls);
if(pAllocd != NULL)
	{
	m_ElSize = ElSize;
	m_NumEls = NumEls;
	m_pAllocd = pAllocd;
	}
else
	m_LastErr = GetLastError();
return(pAllocd);
}

void *
CMemAlloc::Realloc(int NumEls,DWORD HeapFlags)
{
void *pAllocd;
if(m_hHeap == NULL)
	return(NULL);
m_LastErr = 0;
if(NumEls < m_NumEls)		// if actually shrinking then only do so if substantial gains to be made
	{
	if((NumEls * m_ElSize) < 1024)
		NumEls = (1024/m_ElSize) + 1;
	else
		if(NumEls > (m_NumEls*3)/4)
			return(m_pAllocd);
	}
pAllocd = HeapReAlloc(m_hHeap,HeapFlags,m_pAllocd,m_ElSize*NumEls);
if(pAllocd != NULL)
	{
	m_NumEls = NumEls;
	m_pAllocd = pAllocd;
	}
else
	m_LastErr = GetLastError();
return(pAllocd);
}


