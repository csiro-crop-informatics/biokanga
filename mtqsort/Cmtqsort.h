#pragma once

const int cMaxSortThreads = 32;		// allow for a max of this many sort threads

const int cMaxPartStack = 100;		// qsort at most should only require 1 + log2(ElsToSort) stack entries so allow for full 2^64 entries plus a few spare
const int cMergeSortThres = 16;		// switch from qsort to insert sort if <= this number of els to sort

typedef int (*comparer)(const void *, const void *);

#pragma pack(1)

// to avoid recursive calls, unsorted sub-partitions are stacked with effective push/pop through a stack ptr
typedef struct TAG_sSubPartStackEl {
	UINT8 *pCurLeft;						// left boundary for this stacked sub-partition 
	UINT8 *pCurRight;						// right boundary for this stacked sub-partition
} tsSubPartStackEl;

typedef struct TAG_CMTqsort_args {
  class CMTqsort *pThis;
  void *pArray;				// array containing elements to be sorted 
  INT64 NumEls;				// number of elements in array 
  size_t ElSize;			// size in bytes of each element
  comparer CompareFunc;	    // function to compare pairs of elements
} tsCMTqsort_args;

#pragma pack()

class CMTqsort
{
	int m_MaxThreads;							// limit number of threads to be no more than this, defaults to be cMaxSortThreads unless user overrides with call to SetMaxThreads
	int m_CurThreads;							// current number of threads

	tsCMTqsort_args m_ThreadArgs[cMaxSortThreads];	// for holding thread args

	void Exchange(UINT8 *pEl1,			    // exchange this element
	  UINT8 *pEl2,					// with this element
	  size_t ElSize);				// size in bytes of each element

	void InsertSort(UINT8 *pLeft,	// pts to leftmost element		
	    UINT8 *pRight,				// pts to rightmost element
		size_t ElSize,				// size in bytes of each element
		comparer CompareFunc);		// function to compare pairs of elements

	bool							// true if thread was available for handling this partition sort, false if caller needs to do the sort
		ThreadQSort(void *pArray,  	// array containing elements to be sorted
				size_t NumEls,					// number of elements in array
				size_t ElSize,					// size in bytes of each element
				comparer CompareFunc);			// function to compare pairs of elements

	void _mtqsort(bool bUseThreads,						// if true then can create threads to handle sub-partitions
				void *pArray,				// array containing elements to be sorted
				INT64 NumEls,					// number of elements in array
				size_t ElSize,					// size in bytes of each element
				comparer CompareFunc);			// function to compare pairs of elements

#ifdef WIN32
	SRWLOCK m_hRwLock;
	static unsigned int __stdcall _qsort_start (void *args);
#else
	pthread_rwlock_t m_hRwLock;
	static void * _qsort_start (void *args);
#endif

	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);

public:
	CMTqsort(void);
	~CMTqsort(void);

	void SetMaxThreads(int MaxThreads);

	void qsort(void *pArray,					// array containing elements to be sorted
				INT64 NumEls,					// number of elements in array
				size_t ElSize,					// size in bytes of each element
				comparer CompareFunc);			// function to compare pairs of elements
};

