// mtqsort.cpp : Defines the entry point for the console application.
//

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

TRandomCombined<TRanrotWGenerator,TRandomMersenne> RGseeds((int)time(0));

CMTqsort NewSort;

static int
MyStrCmp(const void *pSeq1, const void *pSeq2)
{
return(strcmp((char *)pSeq1, (char *)pSeq2));
}


static int
MyIntCmp(const void *pInt1, const void *pInt2)
{
return(*(int *)pInt1 - *(int *)pInt2);
}

static int
MyInt32Cmp(const void *pInt1, const void *pInt2)
{
return(*(INT32 *)pInt1 - *(INT32 *)pInt2);
}

static int
MyInt64Cmp(const void *pInt1, const void *pInt2)
{
if((*(INT64 *)pInt1 - *(INT64 *)pInt2) < 0)
	return(-1);
else
	if((*(INT64 *)pInt1 - *(INT64 *)pInt2) > 0)
	return(1);
return(0);
}

void printsort (char *s)
{
  char *d = (char*)alloca(strlen(s)+1);
  strcpy(d,s);
  NewSort.qsort(d, (INT64)strlen(d), (size_t)1, MyStrCmp);
  printf("%s -> %s\n", s, d);
}

void test_string ()
{
  printsort("foobarbeque");
  printsort("abcdef");
  printsort("fedcba");
  printsort("abc");
  printsort("acb");
  printsort("bac");
  printsort("bca");
  printsort("cab");
  printsort("cba");
  printsort("ab");
  printsort("ba");
  printsort("a");
  printsort("");
}

void test_int32 ()
{
  UINT32 arr[12], i;
  for (i=0; i<12; i++) arr[i] = 10-i;
  NewSort.qsort(arr, 12, sizeof(UINT32), MyInt32Cmp);
  for (i=0; i<12; i++)
    printf("%d ", arr[i]);
  printf("\n");
}

const int cint64els = 62;

void test_int64 ()
{
UINT64 arr[cint64els], i;

printf("\n In:");
for (i=0; i<cint64els; i++)
	{
	if(i < cint64els*2/3 || i > cint64els-5)
		arr[i] = (rand() * 1001) % 1234567;
	else
		arr[i] = 0;
	printf("%d ", arr[i]);
	}

printf("\nOut:");

NewSort.qsort(arr, cint64els, sizeof(INT64), MyInt64Cmp);
for (i=0; i<cint64els; i++)
  if(i < (cint64els-1) && arr[i] > arr[i+1])
	printf("%d*",arr[i]);
  else
	  printf("%d ", arr[i]);
printf("\n");
}

void benchmark_int64 ()
{
  srand (0);
  UINT64 *arr, i;
  arr = (UINT64*)malloc(10000000*sizeof(UINT64));
  for (i=0; i<10000000; i++) arr[i] = rand();
	NewSort.qsort(arr, 10000000, sizeof(UINT64), MyInt64Cmp);
  free(arr);
}


const size_t cNumInt64s = 800000001; 
const size_t cNumInt32s = 800000000;

void mt_benchmark_int64 ()
{
  srand (0);
  UINT64 *p64;
  UINT64 *arr, Idx;
   printf("\nAllocating mt_benchmark_int64...");
  if((arr = (UINT64*)malloc(cNumInt64s*sizeof(UINT64)))==NULL)
	{
	printf("\nMem allocation failed...");
	exit(1);
	}
  printf("\nSeeding...");

  for (Idx=0; Idx<cNumInt64s; Idx++)
	{
#ifndef USETHISCODE
	if(!(Idx % 100))
		arr[Idx] = Idx;
	else
		if(!(Idx % 33))
			arr[Idx] = 0;
		else
			if(!(Idx % 21))
				arr[Idx] = 100;
			else
				arr[Idx] =  cNumInt64s - Idx;
#else
		arr[Idx] = cNumInt64s - Idx;
#endif
//	  arr[Idx] = 10000;
//	  arr[Idx] = cNumInt64s - Idx;
  //arr[Idx] = RGseeds.IRandom(0,cNumInt64s);
	}

  printf("\nSorting...");

  NewSort.qsort(arr, cNumInt64s, sizeof(UINT64), MyInt64Cmp);
  printf("\nChecking...");

  // now check if have been sorted ascending
  p64 = arr;
  for (Idx=0; Idx<(cNumInt64s-1); Idx++,p64++)
	{
	if(*p64 > p64[1])
		{
		printf("\nNot sorted, Idx %lld == %lld, %lld == %lld ",Idx,*p64,Idx+1,p64[1]);
		break;
		}
	}
  printf("\nChecking mt_benchmark_int64 completed");
  free(arr);
}

void mt_benchmark_int32 ()
{
  srand (0);
  UINT32 *p32;
  UINT32 *arr, Idx;
  printf("\nAllocating for mt_benchmark_int32...");
  if((arr = (UINT32*)malloc(cNumInt32s*sizeof(UINT32)))==NULL)
	  	{
	printf("\nMem allocation failed...");
	exit(1);
	}


  printf("\nSeeding...");
  for (Idx=0; Idx<cNumInt32s; Idx++) 
	  arr[Idx] = 10000;
//	  arr[Idx] = cNumInt32s-Idx;
//	  arr[Idx] = RGseeds.IRandom(0,cNumInt32s);

  printf("\nSorting...");
  NewSort.qsort(arr, (INT64)cNumInt32s, sizeof(UINT32), MyInt32Cmp);
  printf("\nChecking...");

  // now check if have been sorted ascending
  p32 = arr;
  for (Idx=0; Idx<(cNumInt32s-1); Idx++,p32++)
	{
	if(*p32 > p32[1])
		{
		printf("\nNot sorted, Idx %lld == %lld, %lld == %lld ",Idx,*p32,Idx+1,p32[1]);
		break;
		}
	}
  printf("\nChecking mt_benchmark_int32 completed");
  free(arr);
}

int _tmain(int argc, _TCHAR* argv[])
{
  test_string ();
  test_int32 ();
  test_int64 ();
  benchmark_int64 ();
 mt_benchmark_int64();
mt_benchmark_int32();
return 0;
}


















