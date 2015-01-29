#pragma once
// Copyright 2013 CSIRO  ( http://www.csiro.au/ ) 
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License
//   Please contact stuart.stephen@csiro.au for support or 
//   to submit modifications to this source

#include "./pacbiocommon.h"
#include "./PacBioUtility.h"

const int cMinSMRTBellExacts = 20;					// user specified min mumber of exactly matching bases in putative SMRTBell adaptor
const int cDfltMinSMRTBellExacts = 30;				// default mumber of exactly matching bases in putative SMRTBell adaptor
const int cMaxSMRTBellExacts = 46;				    // user specified max mumber of exactly matching bases in putative SMRTBell adaptor

const int cDfltTrim5 = 250;							// default is to after any filtering then to trim 5' end by this many bp
const int cDfltTrim3 = 250;							// default is to after any filtering then to trim 5' end by this many bp
const int cDfltMinReadLen = 5000;					// default is to only report reads of at least this minimum read length after any end trimming

const int cDfltMaxPacBioSeqLen = 0x1ffff;			// default is to allow for PacBio reads of <= 128Kbp
const int cAllocOutBuffSize = 0x0fffff;				// m_pOutBuff allocate output buffer size 

const int cMaxInfiles = 50;							// can accept at most this many input filespecs

typedef enum TAG_ePBPMode {
	ePBPMFilter										// filter and optionally trim PacBio reads
	} etPBPMode;

#pragma pack(1)

typedef struct TAG_sThreadPBFilter {
	int ThreadIdx;					// uniquely identifies this thread
	void *pThis;					// will be initialised to pt to class instance
#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	UINT32 threadID;				// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	int Rslt;						// processing result

	CSSW *pSW;						// Smith-Waterman class
	UINT32 AllocdTargSeqSize;		// current allocation size for buffered target sequence in pTargSeq 	
	etSeqBase *pTargSeq;			// allocated to hold the current probe sequence

	UINT32 AlignErrMem;				// number of times alignments failed because of memory allocation errors
	UINT32 AlignExcessLen;			// number of times alignments failed because length of probe * target was excessive

} tsThreadPBFilter;


#pragma pack()


class CPBFilter
{
	etPBPMode m_PMode;						// processing mode
	int m_MinSMRTBellExacts;				// putative SMRTBell adaptors must contain at least this many exactly matching bases
	int m_SMRTBellFlankSeqLen;				// processing flanking sequences of this length around putative SMRTBell adaptors  
	int m_MinRevCplExacts;					// flanking 5' and RevCpl 3' sequences around putative SMRTBell hairpins must contain at least this many exactly matching bases

	int m_Trim5;							// 5' trim accepted reads by this many bp
	int m_Trim3;							// 3' trim accepted reads by this many bp
	int m_MinReadLen;						// read sequences must be at least this length after any end timming
	int m_NumInputFiles;					// number of input files
	char **m_ppszInputFiles;				// names (wildcards allowed) of input files containing reads to be filtered

	char m_szOutFile[_MAX_PATH];			// name of file in which to write filter accepted and trimmed sequences
	int m_hOutFile;							// handle for file into which to write filter accepted and trimmed sequences
	char *m_pOutBuff;						// used to buffer sequences passing filtering ready to be written to file
	int m_AllocOutBuffSize;					// m_pOutBuff allocated to hold at most this many chars
	int m_OutBuffIdx;						// current number of chars in m_pOutbuff

	int m_NumThreads;						// maximum number of worker threads to use
	int m_TotProcessed;						// total reads processed
	int	m_TotAccepted;						// after filtering accepted this number of reads
	int m_TotPutativeSMRTBells;             // number of putative SMRTBell adaptors
	int	m_TotRejected;                      // rejected this number of reads because of retained PacBio SMRTBell adaptors etc
	int m_TotUnderLen;						// this number reads not accepted because they were underlength

	void Init(void);							// initialise state to that immediately following construction
	void Reset(void);							// reset state
	int LoadTargetSeqs(char *pszTargFile);		// loading sequences form this file 

	int ProcessFiltering(int MaxSeqLen,			// max length sequence expected
							int NumOvlpThreads);	// filtering using at most this many threads

	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);
	void AcquireSerialise(void);
    void ReleaseSerialise(void);
    void AcquireSerialiseMH(void);
	void ReleaseSerialiseMH(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);


#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	HANDLE m_hMtxMHReads;
	SRWLOCK m_hRwLock;
	HANDLE m_hThreadLoadQuerySeqs;
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_mutex_t m_hMtxMHReads;
	pthread_rwlock_t m_hRwLock;
#endif

	static int SortSSWCells(const void *arg1, const void *arg2); // sort putative SMRTBell adaptors by StartTOfs ascending with NumExacts descending
	CPacBioUtility m_PacBioUtility;

public:
	CPBFilter();
	~CPBFilter();

	int PBFilterReads(tsThreadPBFilter *pThreadPar);

	int
	Process(etPBPMode PMode,	// processing mode
		int MinSMRTBellExacts,		// putative SMRTBell adaptors must contain at least this many exactly matching bases
		int SMRTBellFlankSeqLen,    // processing flanking sequences of this length around putative SMRTBell adaptors  
		int MinRevCplExacts,		// flanking 5' and RevCpl 3' sequences around putative SMRTBell hairpins must contain at least this many exactly matching bases
		int Trim5,					// 5' trim accepted reads by this many bp
		int Trim3,					// 3' trim accepted reads by this many bp
		int MinReadLen,				// read sequences must be at least this length after any end timming
		int NumInputFiles,			// number of input files
		char *pszInputFiles[],		// names (wildcards allowed) of input files containing reads to be filtered
		char *pszOutFile,			// name of file in which to write filter accepted and trimmed sequences
		int NumThreads);			// maximum number of worker threads to use

	int
	GenTargIndex(UINT32 MinScaffoldLen,	// individual indexed scaffold sequences must be at least this length
				 char *pszTargFile);	// load sequences from this file

};



