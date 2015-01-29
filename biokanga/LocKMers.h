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

#pragma once

const int cMaxTargCultivarChroms = 50;	// allow the targeted cultivar to have at most this many pseudo chromosomes (each pseudo chrom < 4GB)
const int cMinKMerLen  = 20;			// minimum K-mer length
const int cDfltKMerLen = 50;			// default K-mer length
const int cMaxKMerLen = 100;			// maximum K-mer length
const int cMaxExtKMerLen = 4000;	    // limit length distributions of K-mer lengths out to this maximum - extd K-mers can be much longer than cMaxKMerLen
const int cMinHamming = 1;				// minimum required Hamming separation
const int cDfltHamming = 2;				// default Hamming separation
const int cMaxHamming = 5;				// maximum required Hamming separation

const int cMaxWorkerThreads = 64;		// allow upto 64 threads

const int cBlockReqSize = 0x07ffff;		// each worker thread will request a block containing no more than this total concatenated sequences length
const int cConcatSeqBuffSize = (cBlockReqSize * cMaxWorkerThreads * 8);	// will buffer up to this sized concatenated sequences from which blocks will be allocated

const int cMarkerSeqBuffSize = 0x0ffffff;	// marker sequence buffering used to hold markers prior to writing out to file

// processing modes
typedef enum TAG_ePMode {
	ePMExtdKMers,				// default processing mode
	ePMNoExtdKMers,				// do not attempt to extend K-mers
	ePMPrefixKMers,				// K-mers to share prefix sequence with other cultivars
	ePMplaceholder				// used to set the enumeration range
	} etPMode;

#pragma pack(1)

typedef struct TAG_tsPartialCultivar {
	int PartialCultivarID;		// uniquely identifies this instance (1..N)
	int EntryID;				// suffix array pseudo-chrom identifier
	UINT32 EntryLen;			// pseudo-chrom length
	char szEntryName[cMaxDatasetSpeciesChrom+1]; // pseudo-chrom name
} tsPartialCultivar;

typedef struct TAG_sCultivar {
	tsPartialCultivar Cultivar;		// cultivar details
	UINT32 Status;					// status of this cultivar
} tsCultivar;


typedef struct TAG_sKMerThreadPars {
	int ThreadIdx;						// uniquely identifies this thread
	#ifdef _WIN32
	HANDLE threadHandle;			// handle as returned by _beginthreadex()
	unsigned int threadID;			// identifier as set by _beginthreadex()
#else
	int threadRslt;					// result as returned by pthread_create ()
	pthread_t threadID;				// identifier as set by pthread_create ()
#endif
	class CLocKMers *pThis;			// class
	int Rslt;						// returned result code
	int AllocBlockSeqsSize;			// pBlockSeqs allocated to hold this manybases
	UINT8 *pBlockSeqs;				// used by thread to hold block of sequences for processing
} tsKMerThreadPars;
#pragma pack()

class CLocKMers
{
	CSfxArrayV3 *m_pSfxArray;
	char m_szTargetCultivar[cMaxDatasetSpeciesChrom+1];		// this is the name of the targeted cultivar for which unique K-mers are to be returned
	int m_NumPartialCultivars;		//  this many pseudo chromosome names for targeted cultivar
	int m_CurPartialCultivarID;		//  this partial pseudo-chromosome is currently being buffered for K-mer processing
	UINT32 m_CurPartialCultivarOfs;	// load into buffer starting at this offset
	tsPartialCultivar m_PartialCultivars[cMaxTargCultivarChroms];		// names of pseudo-chromosomes in suffix index file which belong to the targeted cultivar
	int m_CultChromIDs[cMaxTargCultivarChroms];	// partial cultivar pseudo-chroms EntryIDs required by suffix array function MatchesOtherChroms()

	int m_NumSfxEntries;			// total number of cultivar chrom entries in suffix array
	tsCultivar m_AllCultivars[cMaxTargCultivarChroms];	// all cultivars, including those targeted, represented in targeted psudeochromosome sfx array


	int m_PMode;					// processing mode - defaults to maximally extend KMers
	int m_KMerLen;					// this length K-mers
	int m_PrefixLen;				// inter-cultivar shared prefix length
	int m_SuffixLen;				// cultivar specific suffix length
	int m_MinWithPrefix;			// minimum number of cultivars required to have the shared prefix
	int m_MinHamming;				// suffix must be at least this Hamming away from any other K-mer in other species
	int m_NumThreads;				// number of worker threads requested

	// access to following is serialised via AcquireLock()
	UINT32 m_NumBlocks;				// current total number of blocks containing multiple eBaseEOS terminated sequences processed by all threads
	INT64 m_NumKMers;				// current total number of KMers processed by all threads
	UINT32 m_NumPutativeKMers;		// current total number of KMers which are unique to target cultivar
	UINT32 m_NumAcceptedKMers;		// current total number KMers which are unique and more than requested Hammings away in any non-target cultivar
	UINT32 m_NumAcceptedExtdKMers;  // current total number length extended sequences in which any individual subsequence meets the KMer requirement of being more than requested Hammings away in any non-target cultivar  
	// end of serialised access ReleaseLock()

	char m_szDataset[cMaxDatasetSpeciesChrom+1];			// SfxArray dataset name
	bool m_bInitSeqBuffering;		// true if no sequences yet loaded, initialise buffering
	bool m_bAllEntrySeqsLoaded;		// true if all cultivar partial chrom sequences have been loaded into m_pBlockSeqBuff
	bool m_bAllBlocksReturned;		// true if all blocks of sequences have been returned for K-mer markers processing
	UINT32 m_NxtBlockSeqBuffOfs;	// return block of sequences starting at this m_pBlockSeqBuff[m_NxtBlockSeqBuffOfs]
	UINT32 m_BlockSeqBuffLen;		// m_pBlockSeqBuff currently  has been loaded with this many bases
	UINT32 m_AllocBlockSeqBuffSize;	// m_pBlockSeqBuff allocated to hold at most this many bases
	UINT8 *m_pBlockSeqBuff;			// allocated to buffer blocks of sequences to be returned by GetBlockSeqs()

	int m_MaxMarkerLen;				// maximum observed K-mer marker length
	int m_MinMarkerLen;				// maximum observed K-mer marker length
	UINT32 *m_pMarkerLenDist;		// allocated to hold extended K-mer marker length distributions 

	char m_szMarkerFile[_MAX_FNAME];	// write identified markers to this multifasta file
	int m_hOutFile;					// marker output file handle

	bool m_bKMerReads;				// optionally, user can request that all targeted species reads containing the identified K-mer markers also be reported 
	char m_szMarkerReadsFile[_MAX_FNAME];	// write reads containing identified markers to this multifasta file
	int m_hOutReadsFile;			// marker containing reads output file handle

	UINT32 m_MarkerID;				// uniquely identifies marker sequences
	UINT32 m_MarkerBuffOfs;			// concatenate next marker fasta sequence at this offset
	UINT32 m_AllocMarkerBuffSize;	// size of allocated marker buffer
	UINT8 *m_pMarkerBuff;			// allocated to buffer reported marker fasta sequences

	int	FillBlockSeqBuff(void);		// maximally fill sequence buffer	
	
	bool	// returns true if suffix array entry identifier is for a pseudo-chromosome which is part of the target cultivar  
		IsTargetCultivar(int EntryID);		// suffix array entry identifier 

	tsPartialCultivar *				// NULL if EntryID not a targeted cultivar
		LocTargetCultivar(int EntryID);	// suffix array entry identifier

	
	int LocateSpeciesUniqueKMers(tsKMerThreadPars *pPars); // locate all targeted species unique K-mers 

	int LocateSharedUniqueKMers(tsKMerThreadPars *pPars);	// locate all unique K-mers of specified length which are common to all cultivars

	UINT32											// returns number of accepted cultivar specific extended K-mers
		GetKMerProcProgress(UINT32 *pNumBlocks,		// number of blocks processed
					INT64 *pNumKMers,				// these blocks contain this number of K-mers
					UINT32 *pNumPutativeKMers,		// putatively - before check for Hamming - there are this number of unique K-mers mapping to target cultivar
					UINT32 *pNumAcceptedKMers);		// after Hamming check this number of K-mers have been accepted
	
	int
		ReportMarker(int MarkerLen,			// number of bases in marker sequence
			 etSeqBase *pMarkerSeq);		// marker sequence

	int							// marking reads containing the identified marker
		MarkContainingReads(int MarkKMerLen,						// marker is of this length
							   etSeqBase *pMarkStartKMerBase);		// marker sequence

	int
		ReportMarkerRead(int ReadLen,		// number of bases in read sequence
			 etSeqBase *pReadSeq);			// read sequence

#ifdef WIN32
	SRWLOCK m_hRwLock;
	static unsigned int __stdcall KMerThreadStart(void *args);
#else
	pthread_rwlock_t m_hRwLock;
	static void * KMerThreadStart(void *args);
#endif

	void AcquireLock(bool bExclusive);
	void ReleaseLock(bool bExclusive);

	static int MarkersCallback(void *pThis,UINT32 EntryID,etSeqBase *pSeq);
	static int sortpartcultnames(const void *pEl1,const void *pEl2);		// used when sorting cultivar partial pseudo-chrom names

public:
	CLocKMers(void);
	~CLocKMers(void);

	void Reset(bool bSync = false);

	int										// returned block of concatenated sequences total length
		GetBlockSeqs(int MaxLength,			//  maximum total block length to return
						UINT8 *pBlockSeqs);	// copy block of sequences into this buffer


	int
		LocKMers(etPMode PMode,				// processing mode - defaults to 0
		  int KMerLen,					// this length K-mers
	  	  int PrefixLen,				// inter-cultivar shared prefix length
		  int SuffixLen,				// cultivar specific suffix length
		  int MinWithPrefix,			// minimum number of cultivars required to have the shared prefix
		  int MinHamming,				// must be at least this Hamming away from any other K-mer in other cultivars
		  char *pszCultivarName,		// targeted cultivar name for which K-mer markers are required
		  int NumPartialCultivars,		// there are this many pseudo chromosomes for targeted cultivar for which K-mer markers are required
		  char *ppszPartialCultivars[],	// pseudo chromosome names which identify targeted cultivar
		  char *pszSfxPseudoGenome,		// contains pregenerated suffix over psuedochromosomes for each cultivar
		  char *pszMarkerFile,			// output potential markers to this file
		  char *pszMarkerReadsFile,		// optionally output reads containing potential markers to this file
		  int NumThreads);				// max number of threads allowed
};


