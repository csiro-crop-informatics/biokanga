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

const int cMinSMRTBellExacts = 20;					// user specified min mumber of exactly matching bases in putative SMRTBell adapter
const int cDfltMinSMRTBellExacts = 30;				// default mumber of exactly matching bases in putative SMRTBell adapter
const int cMaxSMRTBellExacts = 46;				    // user specified max mumber of exactly matching bases in putative SMRTBell adapter

const int cDfltTrim5 = 100;						// default is to after any filtering then to trim 5' end by this many bp
const int cDfltTrim3 = 100;						// default is to after any filtering then to trim 5' end by this many bp
const int cDfltMinReadLen = 2500;				// default is to only report reads of at least this minimum read length after any end trimming

const int cDfltMaxPacBioSeqLen = 0x1ffff;			// default is to allow for PacBio reads of <= 128Kbp
const int cAllocOutBuffSize = 0x0fffff;				// m_pOutBuff allocate output buffer size 
const int cMaxInfiles = 50;							// can accept at most this many input filespecs

const int cMinContamSeqLen    = 1000;			    // minimum contaminate sequence length accepted - contaminates normally expected to be clone vectors in the 4-8Kbp size range  
const int cMinContamOvlpLen   = 500;			    // requiring contaminants to end overlap by at least this many bp or else be fully contained or containing  
const int cMaxSeedCoreDepth = 10000;					// explore matching cores to this max depth
const int cDeltaCoreOfs = 10;						// offset core windows of coreSeqLen along the probe sequence when checking for overlaps 
const int cCoreSeqLen = 15;							// putative overlaps are explored if there are cores of at least this length in any putative overlap
const int cMinNumCores = 20;						// and if the putative overlap contains at least this many cores
const int cMinPropBinned = 90;						// and if the putative overlap contains at least this proportion (1..100) of 250bp bins binned cores
const int cMaxAcceptHitsPerSeedCore = 500;			// limit accepted hits per seed core to no more this many
const int cDfltMaxProbeSeqLen = cDfltMaxPacBioSeqLen;	// initially allocate for this length probe sequence to be aligned, will be realloc'd as may be required


typedef enum TAG_ePBPMode {
	ePBPMFilter,										// filter and optionally trim PacBio reads
	ePBPMContam                                         // also remove or trim reads containing contaminate sequences
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

	UINT32 SWAlignInstance;			// if contaminate filtering then will be set to the instance of SWAlign to use when aligning to the contaminate sequences
	CSSW *pSW;						// Smith-Waterman class

	UINT32 AlignErrMem;				// number of times alignments failed because of memory allocation errors
	UINT32 AlignExcessLen;			// number of times alignments failed because length of probe * target was excessive

} tsThreadPBFilter;


#pragma pack()


class CPBFilter
{
	etPBPMode m_PMode;						// processing mode
	int m_MinSMRTBellExacts;				// putative SMRTBell adapters must contain at least this many exactly matching bases
	int m_SMRTBellFlankSeqLen;				// processing flanking sequences of this length around putative SMRTBell adapters  
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

	char m_szContamFile[_MAX_PATH];			// name of file to load contaminate sequences from
	int m_ContamARate;						// PacBio expected accuracy event rate, used when contaminate processing
	int m_ContamOvlpLen;					// minimum contaminate overlap length to check for
	CSWAlign *m_pSWAlign;					// used when aligning contaminate sequences

	int m_NumThreads;						// maximum number of worker threads to use
	int m_TotProcessed;						// total reads processed
	int	m_TotAccepted;						// after filtering accepted this number of reads
	int m_TotContamTrimmed;					// this number of accepted reads were contaminate trimmed
	int m_TotPutativeSMRTBells;             // number of putative SMRTBell adapters
	int	m_TotRejected;                      // rejected this number of reads because of retained PacBio SMRTBell adapters
	int m_TotContamRejected;                // rejected this number of reads because of BAC vectors - contained, containing, or overlapping
	int m_TotUnderLen;						// this number reads not accepted because they were underlength

	void Init(void);							// initialise state to that immediately following construction
	void Reset(void);							// reset state

	int ProcessFiltering(int MaxSeqLen,			// max length sequence expected
							int NumOvlpThreads);	// filtering using at most this many threads

	bool m_bMutexesCreated;			// will be set true if synchronisation mutexes have been created
	int CreateMutexes(void);
	void DeleteMutexes(void);
	void AcquireSerialise(void);
    void ReleaseSerialise(void);
	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);


#ifdef _WIN32
	HANDLE m_hMtxIterReads;
	SRWLOCK m_hRwLock;
	HANDLE m_hThreadLoadQuerySeqs;
#else
	pthread_mutex_t m_hMtxIterReads;
	pthread_rwlock_t m_hRwLock;
#endif

	static int SortSSWCells(const void *arg1, const void *arg2); // sort putative SMRTBell adapters by StartTOfs ascending with NumExacts descending
	CPacBioUtility m_PacBioUtility;

public:
	CPBFilter();
	~CPBFilter();

	int PBFilterReads(tsThreadPBFilter *pThreadPar);

	int
	Process(etPBPMode PMode,	// processing mode
		int MinSMRTBellExacts,		// putative SMRTBell adapters must contain at least this many exactly matching bases
		int SMRTBellFlankSeqLen,    // processing flanking sequences of this length around putative SMRTBell adapters  
		int MinRevCplExacts,		// flanking 5' and RevCpl 3' sequences around putative SMRTBell hairpins must contain at least this many exactly matching bases
		int Trim5,					// 5' trim accepted reads by this many bp
		int Trim3,					// 3' trim accepted reads by this many bp
		int MinReadLen,				// read sequences must be at least this length after any end trimming
		int ContamARate,			// PacBio expected accuracy event rate, used when contaminate processing
		int ContamOvlpLen,			// minimum contaminate overlap length to check for
		char *pszContamFile,		// name of file containing contaminate sequences
		int NumInputFiles,			// number of input files
		char *pszInputFiles[],		// names (wildcards allowed) of input files containing reads to be filtered
		char *pszOutFile,			// name of file in which to write filter accepted and trimmed sequences
		int NumThreads);			// maximum number of worker threads to use

};



