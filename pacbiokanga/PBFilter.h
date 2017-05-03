#pragma once
/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "./pacbiocommon.h"
#include "./PacBioUtility.h"

const int cMinSMRTBellExacts = 20;					// user specified min mumber of exactly matching bases in putative SMRTBell adapter
const int cDfltMinSMRTBellExacts = 30;				// default mumber of exactly matching bases in putative SMRTBell adapter
const int cMaxSMRTBellExacts = 46;				    // user specified max mumber of exactly matching bases in putative SMRTBell adapter

const int cDfltTrim5 = 0;						// default is to after any filtering then to trim 5' end by this many bp
const int cDfltTrim3 = 0;						// default is to after any filtering then to trim 5' end by this many bp
const int cDfltMinReadLen = 1000;				// default is to only report reads of at least this minimum read length after any end trimming

const int cDfltMaxPacBioSeqLen = 0x1ffff;			// default is to allow for PacBio reads of <= 128Kbp
const int cAllocOutBuffSize = 0x0fffff;				// m_pOutBuff allocate output buffer size 
const int cMaxInfiles = 200;							// can accept at most this many input filespecs

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
	ePBPMFilter,										// filter and split reads at PacBio SMRTBells and optionally trim PacBio reads
	ePBPMContam                                         // remove or trim reads containing contaminate sequences
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
	int SMRTBellSensitivity; // sensitivity of SMRTBell detection - 5: max, 1: min

	UINT32 AlignErrMem;				// number of times alignments failed because of memory allocation errors
	UINT32 AlignExcessLen;			// number of times alignments failed because length of probe * target was excessive

} tsThreadPBFilter;


#pragma pack()


class CPBFilter
{
	etPBPMode m_PMode;						// processing mode
	int m_Trim5;							// 5' trim accepted reads by this many bp
	int m_Trim3;							// 3' trim accepted reads by this many bp
	int m_MinReadLen;						// read sequences must be at least this length after any end timming
	int m_NumInputFiles;					// number of input files
	char **m_ppszInputFiles;				// names (wildcards allowed) of input files containing reads to be filtered

	char m_szOutFile[_MAX_PATH];			// name of file in which to write filter accepted and trimmed sequences
	int m_hOutFile;							// handle for file into which to write filter accepted and trimmed sequences
	char *m_pszOutBuff;						// used to buffer sequences passing filtering ready to be written to file
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
	int m_TotSMRTBells;						// number of SMRTBell adapters located
	int m_TotSMRTBellSubSeqs;				// sequences with SMRTBell adapters trimmed into this many subsequences
	int	m_TotRejected;                      // rejected this number of reads because of retained PacBio SMRTBell adapters
	int m_TotContamRejected;                // rejected this number of reads because of BAC vectors - contained, containing, or overlapping
	int m_TotUnderLen;						// this number reads not accepted because they were under length

	void Init(void);							// initialise state to that immediately following construction
	void Reset(void);							// reset state

	int ProcessFiltering(int SMRTBellSensitivity, // sensitivity of SMRTBell detection - 5: max, 1: min
						int MaxSeqLen,			// max length sequence expected
						int NumOvlpThreads);	// filtering using at most this many threads

	int			// eBSFSuccess or error code if write failed
		WriteFastaFile(int SeqID,	// primary sequence identifier
			int SecSeqID,			// secondary sequence identifier
			int SeqLen,				// sequence to write is this length
			etSeqBase *pSeq);		// sequence

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
		int SMRTBellSensitivity, // sensitivity of SMRTBell detection - 5: max, 1: min
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



