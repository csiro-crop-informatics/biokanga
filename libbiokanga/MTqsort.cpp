/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "../libbiokanga/commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "../libbiokanga/commhdrs.h"
#endif

// need for speed rather than space...
#pragma optimize("t", on)

// constructor
CMTqsort::CMTqsort(void)
{
m_MaxThreads = cDfltSortThreads;
m_CurThreads = 0;
memset(m_ThreadArgs,0,sizeof(m_ThreadArgs));
#ifdef _WIN32
InitializeSRWLock(&m_hRwLock);
#else
pthread_rwlock_init (&m_hRwLock,NULL);
#endif
}

// destructor
CMTqsort::~CMTqsort(void)
{
#ifndef _WIN32
pthread_rwlock_destroy(&m_hRwLock);
#endif
}

// AcquireLock
// Aquire lock, either exclusive or read only, on class shared instance vars
void 
CMTqsort::AcquireLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	AcquireSRWLockExclusive(&m_hRwLock);
else
	AcquireSRWLockShared(&m_hRwLock);
#else
if(bExclusive)
	pthread_rwlock_wrlock(&m_hRwLock);
else
	pthread_rwlock_rdlock(&m_hRwLock);
#endif
}


// ReleaseLock
// Release lock, must be exclusive or read only as used when AcquireLock was used
void
CMTqsort::ReleaseLock(bool bExclusive)
{
#ifdef _WIN32
if(bExclusive)
	ReleaseSRWLockExclusive(&m_hRwLock);
else
	ReleaseSRWLockShared(&m_hRwLock);
#else
pthread_rwlock_unlock(&m_hRwLock);
#endif
}

// SetMaxThreads
// Sets maximum number of threads to use, if 0 then resets to cMaxSortThreads
void 
CMTqsort::SetMaxThreads(int MaxThreads)
{
if(MaxThreads <= 0 || MaxThreads > cMaxSortThreads)
	MaxThreads = cMaxSortThreads;
AcquireLock(true);
m_MaxThreads = MaxThreads;
ReleaseLock(true);
}


// _qsort_start
// Thread start - simply unpacks it's args into a call to mtqsort
#ifdef WIN32
unsigned int __stdcall CMTqsort::_qsort_start (void *args)
{
#else
void * CMTqsort::_qsort_start (void *args)
{
pthread_detach(pthread_self());
#endif
tsCMTqsort_args *pArgs = (tsCMTqsort_args *)args;

// use library qsort as should be more optimised than this classes _mtqsort so use it..
if(pArgs->NumEls < (cMinUseLibQsort * 2))
	::qsort(pArgs->pArray, (size_t)pArgs->NumEls, pArgs->ElSize, pArgs->CompareFunc);
else
	pArgs->pThis->_mtqsort(true,pArgs->pArray, pArgs->NumEls, pArgs->ElSize, pArgs->CompareFunc);

// one less thread sorting, allow other threads to start up
pArgs->pThis->AcquireLock(true);
if(pArgs->pThis->m_CurThreads > 0)
	pArgs->pThis->m_CurThreads -= 1;
pArgs->ElSize = 0;						// releases for other threads to use
pArgs->pThis->ReleaseLock(true);
#ifdef WIN32
ExitThread(1);
#else
return NULL;
#endif
}



// Exchange
// Exchange two elements each of size ElSize, note it is assumed that elements are not overlapped!
void 
CMTqsort::Exchange(UINT8 *pEl1,			    // exchange this element
	  UINT8 *pEl2,					// with this element
	  size_t ElSize)				// size in bytes of each element
{
size_t Idx;

union {
	UINT8 T8;
	UINT16 T16;
	UINT32 T32;
	UINT64 T64;
} Tmp;

if(pEl1 == pEl2 || ElSize < 1)
	return;

switch(ElSize) {
	case 1:
		Tmp.T8 = *pEl1;
		*pEl1 = *pEl2;
		*pEl2 = Tmp.T8;
		return;
	case 2:
		Tmp.T16 = *(UINT16 *)pEl1;
		*(UINT16 *)pEl1 = *(UINT16 *)pEl2;
		*(UINT16 *)pEl2 = Tmp.T16;
		return;

	case 3:
		Tmp.T16 = *(UINT16 *)pEl1;
		*(UINT16 *)pEl1 = *(UINT16 *)pEl2;
		*(UINT16 *)pEl2 = Tmp.T16;
		pEl1 += sizeof(UINT16);
		pEl2 += sizeof(UINT16);
		Tmp.T8 = *pEl1;
		*pEl1 = *pEl2;
		*pEl2 = Tmp.T8;
		return;

	case 4:
		Tmp.T32 = *(UINT32 *)pEl1;
		*(UINT32 *)pEl1 = *(UINT32 *)pEl2;
		*(UINT32 *)pEl2 = Tmp.T32;
		return;

	case 5:
		Tmp.T32 = *(UINT32 *)pEl1;
		*(UINT32 *)pEl1 = *(UINT32 *)pEl2;
		*(UINT32 *)pEl2 = Tmp.T32;
		pEl1 += sizeof(UINT32);
		pEl2 += sizeof(UINT32);
		Tmp.T8 = *pEl1;
		*pEl1 = *pEl2;
		*pEl2 = Tmp.T8;
		return;

	case 6:
		Tmp.T32 = *(UINT32 *)pEl1;
		*(UINT32 *)pEl1 = *(UINT32 *)pEl2;
		*(UINT32 *)pEl2 = Tmp.T32;
		pEl1 += sizeof(UINT32);
		pEl2 += sizeof(UINT32);
		Tmp.T16 = *(UINT16 *)pEl1;
		*(UINT16 *)pEl1 = *(UINT16 *)pEl2;
		*(UINT16 *)pEl2 = Tmp.T16;
		return;

	case 7:
		Tmp.T32 = *(UINT32 *)pEl1;
		*(UINT32 *)pEl1 = *(UINT32 *)pEl2;
		*(UINT32 *)pEl2 = Tmp.T32;
		pEl1 += sizeof(UINT32);
		pEl2 += sizeof(UINT32);
		Tmp.T16 = *(UINT16 *)pEl1;
		*(UINT16 *)pEl1 = *(UINT16 *)pEl2;
		*(UINT16 *)pEl2 = Tmp.T16;
		pEl1 += sizeof(UINT16);
		pEl2 += sizeof(UINT16);
		Tmp.T8 = *pEl1;
		*pEl1 = *pEl2;
		*pEl2 = Tmp.T8;
		return;

	case 8:
		Tmp.T64 = *(UINT64 *)pEl1;
		*(UINT64 *)pEl1 = *(UINT64 *)pEl2;
		*(UINT64 *)pEl2 = Tmp.T64;
		return;

	default:					// some strange size, first exchange as UINT32's then remainder as UINT8s
		for (Idx=0; Idx < (UINT32)(ElSize - ElSize % sizeof(UINT32)); Idx += sizeof(UINT32)) 
			{
			Tmp.T32 = *(UINT32 *)pEl1;
			*(UINT32 *)pEl1 = *(UINT32 *)pEl2;
			*(UINT32 *)pEl2 = Tmp.T32;
			pEl1 += sizeof(UINT32);
			pEl2 += sizeof(UINT32);
			}

		  for (; Idx < ElSize; Idx++) 
			{
			Tmp.T8 = *pEl1;
			*pEl1++ = *pEl2;
			*pEl2++ = Tmp.T8;
			}
		return;
	}
}

// InsertSort
// Insertion sort used when the number of elements in partition is below cInsertSortMinLen
void
CMTqsort::InsertSort(UINT8 *pLeft,	// pts to leftmost element		
	    UINT8 *pRight,				// pts to rightmost element
		size_t ElSize,				// size in bytes of each element
		comparer CompareFunc)		// function to compare pairs of elements
{
UINT8 *pProbe;
UINT8 *pMax;

while (pRight > pLeft) {
	pMax = pLeft;
    for (pProbe = pLeft + ElSize; pProbe <= pRight; pProbe += ElSize)
	   if (CompareFunc( pProbe, pMax) > 0)
                pMax = pProbe;
    Exchange(pMax, pRight, ElSize);
	pRight -= ElSize;
    }
}





bool							// true if thread was available for handling this partition sort, false if caller needs to do the sort
CMTqsort::ThreadQSort(void *pArray,
				size_t NumEls,
				size_t ElSize,
				comparer CompareFunc)
{
int Idx;
tsCMTqsort_args *pArgs;
#ifdef _WIN32
HANDLE threadHandle;			// handle as returned by _beginthreadex()
unsigned int threadID;			// identifier as set by _beginthreadex()
#else
int threadRslt;					// result as returned by pthread_create ()
pthread_t threadID;				// identifier as set by pthread_create ()
#endif

AcquireLock(true);
if(m_CurThreads < m_MaxThreads)
	{
	m_CurThreads += 1;			// will be starting a new thread
	pArgs = m_ThreadArgs;
	for(Idx = 0; Idx < m_MaxThreads; Idx++,pArgs++)
		{
		if(pArgs->ElSize == 0)	// this thread arg is currently avail for reuse
			break;
		}
	if(Idx == m_MaxThreads)		// double check there was a free thread arg 
		{
		ReleaseLock(true);
		return(false);
		}
	pArgs->pThis = this;
	pArgs->NumEls = NumEls;
	pArgs->pArray = pArray;
	pArgs->ElSize = ElSize;
	pArgs->CompareFunc = CompareFunc;
	ReleaseLock(true);

#ifdef WIN32
	threadHandle = (HANDLE)_beginthreadex(NULL,0x0fffff,_qsort_start,(void*)pArgs,0,&threadID);
#else
    threadRslt = pthread_create(&threadID, NULL , _qsort_start, (void*)pArgs);
#endif
	return(true);
	}
ReleaseLock(true);
return(false);
}

void 
CMTqsort::_mtqsort (bool bUseThreads,		// if true then can create threads to handle sub-partitions
				void *pArray,				// array containing elements to be sorted
				INT64 NumEls,					// number of elements in array
				size_t ElSize,					// size in bytes of each element
				comparer CompareFunc)			// function to compare pairs of elements
{
UINT8 *pCurPartStart;							// current partition start
UINT8 *pCurPartEnd;								// current partition end
UINT8 *pCurPartMid;								// pivot point in current partition, hopefully will be the median value so partition will be split into two equal sized sub partitions!
UINT8 *pLow;									// used when traversing partition starting from start towards end 
UINT8 *pHigh;									// used when traversing partition starting from end towards start 
size_t NumElsCurPart;							// number of elements in current sub-partition
tsSubPartStackEl SubPartStackEls[cMaxPartStack];		// stack of sub-partitions to yet to be processed 
tsSubPartStackEl *pSubPartStackEl;				// stack ptr

if (NumEls < 2)									// anything to sort?		
    return;                

pSubPartStackEl = &SubPartStackEls[0];                 
pCurPartStart = (UINT8 *)pArray;
pCurPartEnd = (UINT8 *)pArray + ElSize * (NumEls-1);    

recurse:
	{
    NumElsCurPart = (pCurPartEnd - pCurPartStart) / ElSize + 1;        
	if(!bUseThreads || !((NumEls/NumElsCurPart >= 4) && NumElsCurPart > (cMinUseLibQsort * 2) && 
		NumElsCurPart < cMaxUseThreadQsort &&
		ThreadQSort(pCurPartStart,NumElsCurPart,ElSize, CompareFunc)))
		{
		if (NumElsCurPart <= cMergeSortThres)										// with small number of els then more efficent to do a insert sort than continueing with the qsort
			InsertSort(pCurPartStart, pCurPartEnd, ElSize, CompareFunc);
		else 
			 {		
			 // select a pivot as being the median of 3, with luck this may result in a near even split of the current partition
			 pCurPartMid = pCurPartStart + (NumElsCurPart / 2) * ElSize;
			 if (CompareFunc( pCurPartStart, pCurPartMid) > 0)
				Exchange(pCurPartStart, pCurPartMid, ElSize);
			if (CompareFunc( pCurPartStart, pCurPartEnd) > 0)
				Exchange(pCurPartStart, pCurPartEnd, ElSize);
			if (CompareFunc( pCurPartMid, pCurPartEnd) > 0)
				Exchange(pCurPartMid, pCurPartEnd, ElSize);
        
			pLow = pCurPartStart;
			pHigh = pCurPartEnd;

			while(1) 
				{
				 if (pCurPartMid > pLow) 
					do  {
						pLow += ElSize;
						} 
					while (pLow < pCurPartMid && CompareFunc( pLow, pCurPartMid) <= 0);
				

				if (pCurPartMid <= pLow)
					do  {
						pLow += ElSize;
						} while (pLow <= pCurPartEnd && CompareFunc( pLow, pCurPartMid) <= 0);
				
				do  {
					pHigh -= ElSize;
					} 
				while (pHigh > pCurPartMid && CompareFunc(pHigh, pCurPartMid) > 0);
			
				if (pHigh < pLow)
					break;

				Exchange(pLow, pHigh, ElSize);

				if (pCurPartMid == pHigh)
					pCurPartMid = pLow;
				}


			pHigh += ElSize;
			if (pCurPartMid < pHigh)
				do  {
					pHigh -= ElSize;
					} 
				while (pHigh > pCurPartMid && CompareFunc(pHigh, pCurPartMid) == 0);
			
			if (pCurPartMid >= pHigh)
				do  {
					pHigh -= ElSize;
					} 
				while (pHigh > pCurPartStart && CompareFunc(pHigh, pCurPartMid) == 0);
			

			if (pHigh - pCurPartStart >= pCurPartEnd - pLow ) 
				{
				if (pCurPartStart < pHigh) 
					{
					pSubPartStackEl->pCurLeft = pCurPartStart;
					pSubPartStackEl->pCurRight = pHigh;
					pSubPartStackEl += 1;
 					}                   

				if (pLow < pCurPartEnd) 
					{
					pCurPartStart = pLow;
					goto recurse;       
					}
				}
			else 
				{
				if(pLow < pCurPartEnd) 
					{
					pSubPartStackEl->pCurLeft = pLow;
					pSubPartStackEl->pCurRight = pCurPartEnd;
					pSubPartStackEl += 1;
					}

				if(pCurPartStart < pHigh) 
					{
					pCurPartEnd = pHigh;
					goto recurse;           
					}
				}
			}
		}
	}

pSubPartStackEl -= 1;
if (pSubPartStackEl >= &SubPartStackEls[0]) 
	{
	pCurPartStart = pSubPartStackEl->pCurLeft;
	pCurPartEnd = pSubPartStackEl->pCurRight;
	goto recurse; 
	}

return;       
}

// threaded qsort 
void 
CMTqsort::qsort(void *pArray,
				INT64 NumEls,
				size_t ElSize,
				comparer CompareFunc)
{
int CurThreads;
if(pArray == NULL || NumEls <= 1 || ElSize < 1 || CompareFunc == NULL)
	return;

if(NumEls < cMinUseLibQsort)
	return(::qsort(pArray,(size_t)NumEls,ElSize,CompareFunc));

m_CurThreads = 1;		// this thread counts as the first, additional threads will be created to process sub-partitions up the limit of m_MaxThreads
memset(m_ThreadArgs,0,sizeof(m_ThreadArgs));

// start sort processing
_mtqsort(true, pArray,NumEls,ElSize,CompareFunc);

// spin-wait for all threads which may still be processing sub-partitions to complete
do {
	AcquireLock(false);
	CurThreads = m_CurThreads;
	ReleaseLock(false);
	if(CurThreads > 1)
#ifdef _WIN32
		Sleep(500);
#else
		sleep(1);
#endif	
	}
while(CurThreads > 1);
}













