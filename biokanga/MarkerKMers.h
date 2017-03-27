/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#pragma once

const int cMinKMerLen  = 20;			// minimum K-mer length
const int cDfltKMerLen = 50;			// default K-mer length
const int cMaxKMerLen = 200;			// maximum K-mer length

const int cMaxWorkerThreads = 128;			// limiting max number of threads to this many

const int cMarkerSeqBuffSize = 0x0fffff;	// marker sequence buffering used to hold markers prior to writing out to file
const int cAllocNumPutativeSeqs = 0x0fffff;	// allocate for this many putative marker sequences, and realloc in this many increments as may be required

const UINT8 cMarkerDupFlg = 0x01;			// marker prefix sequence is a duplicate
const UINT8 cMarkerOvlFlg = 0x02;			// marker prefix sequence overlaps onto another prefix sequence
const UINT8 cMarkerAntiFlg = 0x04;			// marker prefix sequence is antisense to another prefix sequence


// processing modes
typedef enum TAG_ePMode {
	ePMSenseAntiKMers,			// default processing mode is for both sense and antisense K-mer processing
	ePMNSenseKMers,				// process sense K-mers only
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

#pragma pack(1)

typedef struct TAG_sCultivar {
	int EntryID;				// suffix array pseudo-chrom (cultivar) identifier
	UINT32 EntryLen;			// pseudo-chrom length
	char szEntryName[cMaxDatasetSpeciesChrom+1]; // pseudo-chrom name
	UINT32 Status;					// status of this cultivar
} tsCultivar;

typedef struct TAG_sPutMarker {
	UINT32 MarkerID;		// uniquely identifies this marker
	UINT16  NumCultivars;	// number of cultivars in which marker was located
	UINT32 SenseCnts;		// number of marker K-Mers on sense strand
	UINT32 AntisenseCnts;	// number of marker K-Mers on antisense strand
	UINT8  Flags;			// set to a combination of cMarkerDupFlg, cMarkerOvlFlg, cMarkerAntiFlg  
	UINT8  MarkerSeq[1];	// to hold marker sequence bases
} tsPutMarker;

typedef struct TAG_sKMerThreadPars {
	int ThreadIdx;						// uniquely identifies this thread
	#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	class CMarkerKMers *pThis;		// class instance
	INT64 StartSfxIdx;				// thread to process from this suffix index inclusive
	INT64 EndSfxIdx;				// thread to process until this suffix index inclusive
	int Rslt;						// returned result code
} tsKMerThreadPars;
#pragma pack()

class CMarkerKMers
{
	CSfxArrayV3 *m_pSfxArray;
	int m_NumSfxEntries;			// total number of cultivar chrom entries in suffix array
	tsCultivar m_AllCultivars[cMaxCultivars];	// all cultivars represented in targeted psudeochromosome sfx array


	int m_PMode;					// processing mode - defaults to proccessing K-Mers sense and antisense
	int m_KMerLen;					// this length K-mers
	int m_PrefixLen;				// K-mer prefix length
	int m_SuffixLen;				// K-mer suffix length
	int m_MinWithPrefix;			// minimum number of cultivars required to have the shared prefix
	int m_MaxHomozygotic;			// only report prefixes if K-Mer suffixes are homozygotic between a maximum of this many cultivars, if 0 then no homozygotic check
	int m_NumThreads;				// number of worker threads requested

	UINT32 m_NumPrefixKMers;		// current total number of accepted unique prefixed KMers
	INT64 m_TotSenseCnts;			// current total number of accepted KMers on sense strand
	INT64 m_TotAntisenseCnts;		// current total number of accepted KMers on antisense strand
	
	char m_szDataset[cMaxDatasetSpeciesChrom+1];			// SfxArray dataset name

	char m_szMarkerFile[_MAX_FNAME];	// write identified markers to this multifasta file
	int m_hOutFile;					// marker output file handle

	UINT32 m_MarkerID;				// uniquely identifies marker sequences

	UINT32 m_MarkerBuffOfs;			// concatenate next marker fasta sequence at this offset
	UINT32 m_AllocMarkerBuffSize;	// size of allocated marker buffer
	UINT8 *m_pMarkerBuff;			// allocated to buffer reported marker fasta sequences

	int m_PutMarkerSize;			// size of each tsPutMarker allowing for m_PrefixLen sequence bases 
	UINT32 m_NumPutMarkers;			// current number of putative marker sequences in m_pPutativeSeqs
	size_t m_AllocPutMarkersSize;	// current memory allocation size of m_pPutativeSeqs
	tsPutMarker *m_pPutMarkers;		// used to hold all putative marker sequences
	size_t m_AllocPutMarkersIndexSize; // current memory allocation size of m_pPutativeSeqs
	UINT32 *m_pPutMarkersIndex;		// used to hold sorted index over m_pPutMarkers

	int LocateSharedPrefixKMers(tsKMerThreadPars *pPars); // locate K-mers with shared prefixes 

	int LocateSharedUniqueKMers(tsKMerThreadPars *pPars);	// locate all unique K-mers of specified length which are common to all cultivars

	bool											// true if overlapping by m_PrefixLen-1 onto at least one other prefix marker
		IsMarkerOverlapping(tsPutMarker *pMarker);	// marker to check if overlapping onto any other sequence

	tsPutMarker *										// NULL if no others antisense, else the other antisense marker
		GetMarkerAntisense(tsPutMarker *pMarker);		// marker to check if antisense to any other marker sequence
							

	INT64											// returns number of K-Mers processed
		GetKMerProcProgress(INT64 *pTotSenseCnts,	// number of K-Mers on sense strand
					INT64 *pTotAntisenseCnts);		// number of K-Mers on antisense strand
	
	int
		ReportMarker(tsPutMarker *pMarker);			// marker to report


#ifdef WIN32
	SRWLOCK m_hRwLock;
	CRITICAL_SECTION m_hSCritSect;
	static unsigned int __stdcall KMerThreadStart(void *args);
#else
	pthread_rwlock_t m_hRwLock;
	pthread_spinlock_t m_hSpinLock;
	static void * KMerThreadStart(void *args);
#endif

	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);

	void EnterCritSect(void);
	void LeaveCritSect(void);

	static int MarkersCallback(void *pThis,tsKMerCultsCnts *pCultsCnts);
	static int sortpartcultpseudonames(const void *pEl1,const void *pEl2);		// used when sorting cultivar partial pseudo-chrom names

	static int SortPutativeSeqs(const void *arg1, const void *arg2);			// used when sorting putative marker sequences
	static int SortNumCultivarsCnts(const void *arg1, const void *arg2);		// used when sorting putative markers by NumCultivars and sense/antisense counts

public:
	CMarkerKMers(void);
	~CMarkerKMers(void);

	void Reset(bool bSync = false);

	int										// returned block of concatenated sequences total length
		GetBlockSeqs(int MaxLength,			//  maximum total block length to return
						UINT8 *pBlockSeqs);	// copy block of sequences into this buffer


	int
		LocKMers(etPMode PMode,			// processing mode - defaults to 0
		  int KMerLen,					// this length K-mers
	  	  int PrefixLen,				// inter-cultivar shared prefix length
		  int SuffixLen,				// cultivar specific suffix length
		  int MinWithPrefix,			// minimum number of cultivars required to have the shared prefix
		  int MaxHomozygotic,			// only report prefixes if K-Mer suffixes are homozygotic between a maximum of this many cultivars, if 0 then no homozygotic check
		  char *pszSfxPseudoGenome,		// contains pregenerated suffix over psuedochromosomes for each cultivar
		  char *pszMarkerFile,			// output potential markers to this file
		  int NumThreads);				// max number of threads allowed
};


