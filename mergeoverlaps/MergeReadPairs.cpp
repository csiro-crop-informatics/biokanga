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

#include "MergeReadPairs.h"


CMergeReadPairs::CMergeReadPairs(void)
{
m_pszMSeqs = NULL;
m_pszUnmergedP1Seqs = NULL;
m_pszUnmergedP2Seqs = NULL;
m_hIn5Reads = -1;
m_hOutMerged = -1;
m_hOut5Unmerged = -1;
m_hOut3Unmerged = -1;

m_ProcPhase = ePPUninit;
Reset(false);
}


CMergeReadPairs::~CMergeReadPairs(void)
{
Reset(false);
}

void
CMergeReadPairs::Reset(bool bSync)
{
if(m_ProcPhase != ePPUninit)
	{
	if(m_hIn5Reads != -1)
		{
		close(m_hIn5Reads);
		m_hIn5Reads = -1;
		}

	if(m_hIn3Reads != -1)
		{
		close(m_hIn3Reads);
		m_hIn3Reads = -1;
		}

	if(m_pszMSeqs != NULL)
		{
		if(bSync && m_CurMSeqLen && m_hOutMerged != -1)
			CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,m_CurMSeqLen);
		delete m_pszMSeqs;
		m_pszMSeqs = NULL;
		}
	m_AllocdMSeqs = 0;

	if(m_pszUnmergedP1Seqs != NULL)
		{
		if(bSync && m_CurUnmergedP1Seqs && m_hOut5Unmerged != -1)
			CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,m_CurUnmergedP1Seqs);
		delete m_pszUnmergedP1Seqs;
		m_pszUnmergedP1Seqs = NULL;
		}
	m_AllocdUnmergedP1Seqs = 0;

	if(bSync && m_pszUnmergedP2Seqs != NULL)
		{
		if(m_CurUnmergedP2Seqs && m_hOut3Unmerged != -1)
			CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,m_CurUnmergedP2Seqs);
		delete m_pszUnmergedP2Seqs;
		m_pszUnmergedP2Seqs = NULL;
		}
	m_AllocdUnmergedP2Seqs = 0;

	if(m_hOutMerged != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOutMerged);
#else
			fsync(m_hOutMerged);
#endif
		close(m_hOutMerged);
		m_hOutMerged = -1;
		}

	if(m_hOut5Unmerged != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOut5Unmerged);
#else
			fsync(m_hOut5Unmerged);
#endif
		close(m_hOut5Unmerged);
		m_hOut5Unmerged = -1;
		}

	if(m_hOut3Unmerged != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOut3Unmerged);
#else
			fsync(m_hOut3Unmerged);
#endif
		close(m_hOut3Unmerged);
		m_hOut3Unmerged = -1;
		}
	}
else
	{
	m_hIn5Reads = -1;			// file handle for opened 5' reads
	m_hIn3Reads = -1;			// file handle for opened 3' reads
	m_hOutMerged = -1;			// file handle for output merged read microcontigs
	m_hOut5Unmerged = -1;		// file handle for output 5' unmerged reads
	m_hOut3Unmerged = -1;		// file handle for output 3' unmerged reads
	m_pszMSeqs = NULL;			// will be allocated for holding merged sequences ready for writing to file
	m_pszUnmergedP1Seqs = NULL;  // will be allocated for holding unmerged P1 sequences ready for writing to file
	m_pszUnmergedP2Seqs = NULL;  // will be allocated for holding unmerged P2 sequences ready for writing to file
	m_AllocdMSeqs = 0;
	m_AllocdUnmergedP1Seqs = 0;
	m_AllocdUnmergedP2Seqs = 0;
	}

m_szIn5ReadsFile[0] = 0;	// 5' reads are in this file
m_szIn3ReadsFile[0] = 0;	// 3' reads are in this file
m_szMergeOutFile[0] = 0;	// write merged overlaping reads to this file
m_sz5UnmergedOutFile[0] = 0;	// write 5' unmerged reads to this file
m_sz3UnmergedOutFile[0] = 0;   // write 3' unmerged reads to this file

m_PMode = ePMdefault;
m_OFormat = eOFauto;		
m_bIsFastq = false;			// if true then input reads are fastq
m_bAppendOut = false;		// if true then append if output files exist, otherwise trunctate existing output files 
m_MinOverlap = 0;			// reads must overlap by at least this many bases
m_MaxOverlapPropSubs = 0;	// and overlap can have at most this proportion of substitutions; if non-zero then a floor of 1 sub is allowed

m_CurMSeqLen = 0;			// m_pszMSeqs currently holds this many chars
m_CurUnmergedP1Seqs = 0;	// m_pszUnmergedP1Seqs currently holds this many chars
m_CurUnmergedP2Seqs = 0;	// m_pszUnmergedP2Seqs currently holds this many chars
m_ProcPhase = ePPReset;
}

// OpenFiles
// Initialise and open files for processing
int
CMergeReadPairs::OpenFiles(void)
{

// check if files specified for open/create/append...
if(m_ProcPhase != ePPRdyOpen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"OpenFiles: unexpected state error");
	return(eBSFerrInternal);
	}

#ifdef _WIN32
m_hIn5Reads = open(m_szIn5ReadsFile, O_READSEQ ); // file access is sequential..
#else
m_hIn5Reads = open64(m_szIn5ReadsFile, O_READSEQ ); 
#endif

if(m_hIn5Reads == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s for 5' PE processing - %s",m_szIn5ReadsFile,strerror(errno));
	Reset(false);
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opened %s for 5' PE processing",m_szIn5ReadsFile);

#ifdef _WIN32
m_hIn3Reads = open(m_szIn3ReadsFile, O_READSEQ ); // file access is sequential..
#else
m_hIn3Reads = open64(m_szIn3ReadsFile, O_READSEQ ); 
#endif

if(m_hIn3Reads == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s for 5' PE processing - %s",m_szIn3ReadsFile,strerror(errno));
	Reset(false);
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Opened %s for 3' PE processing",m_szIn5ReadsFile);


if(m_bAppendOut)
	{
#ifdef _WIN32
	m_hOutMerged = open(m_szMergeOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
	m_hOutMerged = open(m_szMergeOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
	if(m_hOutMerged == -1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_szMergeOutFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	// seek to end of file ready for appending
	lseek(m_hOutMerged,0,SEEK_END);

	if(m_sz5UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
		if(m_hOut5Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz5UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		// seek to end of file ready for appending
		lseek(m_hOut5Unmerged,0,SEEK_END);
		}

	if(m_sz3UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
		if(m_hOut3Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output 3' unmerged file '%s'",m_sz3UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		// seek to end of file ready for appending
		lseek(m_hOut3Unmerged,0,SEEK_END);
		}
	}
else
	{
#ifdef _WIN32
	m_hOutMerged = open(m_szMergeOutFile,O_CREATETRUNC );
#else
	if((m_hOutMerged = open(m_szMergeOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutMerged,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_szMergeOutFile,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif
	if(m_hOutMerged == -1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_szMergeOutFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	lseek(m_hOutMerged,0,SEEK_SET);

	if(m_sz5UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_CREATETRUNC );
#else
		if((m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOut5Unmerged,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_sz5UnmergedOutFile,strerror(errno));
					Reset(false);
					return(eBSFerrCreateFile);
					}
#endif
		if(m_hOut5Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz5UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		lseek(m_hOut5Unmerged,0,SEEK_SET);
		}

	if(m_sz3UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_CREATETRUNC );
#else
		if((m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOut3Unmerged,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_sz3UnmergedOutFile,strerror(errno));
					Reset(false);
					return(eBSFerrCreateFile);
					}
#endif
		if(m_hOut3Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz3UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		lseek(m_hOut3Unmerged,0,SEEK_SET);
		}
	}
m_ProcPhase = ePRdyProc;
return(eBSFSuccess);
}

int
CMergeReadPairs::MergeOverlaps(etPMode PMode,	// processing mode 
			  etOFormat OFormat,				// output file format
			  bool bAppendOut,					// if true then append if output files exist, otherwise trunctate existing output files
			  int StartNum,					// initial starting sequence identifier
			  int MinOverlap,				// reads must overlap by at least this many bases
              int MaxOverlapPropSubs,		// and overlap can have at most this proportion of substitutions, if > 0 then floor of 1 sub allowed
			  char *pszIn5ReadsFile,		// 5' reads are in this file
			  char *pszIn3ReadsFile,		// 3' reads are in this file
			  char *pszMergeOutFile)		// write merged overlaping reads to this file
{
int Rslt;

if(MinOverlap < 1 ||					// must be at least a 1 base overlap required!
   MaxOverlapPropSubs < 0 ||			// must be reasonable number of allowed substitutions in overlap
   pszIn5ReadsFile == NULL || pszIn5ReadsFile[0] == '\0' ||
   pszIn3ReadsFile == NULL || pszIn3ReadsFile[0] == '\0' ||
   pszMergeOutFile == NULL || pszMergeOutFile[0] == '\0')
   	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"MergeOverlaps: inconsistency in function parameters");
	return(eBSFerrParams);
	}
Reset(false);	
m_PMode = PMode;
m_OFormat = OFormat;

m_bAppendOut = bAppendOut;
m_MinOverlap = MinOverlap;
m_MaxOverlapPropSubs = MaxOverlapPropSubs;
strncpy(m_szIn5ReadsFile,pszIn5ReadsFile,_MAX_PATH);
m_szIn5ReadsFile[_MAX_PATH-1] = '\0';
strncpy(m_szIn3ReadsFile,pszIn3ReadsFile,_MAX_PATH);
m_szIn3ReadsFile[_MAX_PATH-1] = '\0';
strncpy(m_szMergeOutFile,pszMergeOutFile,_MAX_PATH);
m_szMergeOutFile[_MAX_PATH-1] = '\0';
if(PMode == ePMseparate)
	{
	strncpy(m_sz5UnmergedOutFile,m_szMergeOutFile,_MAX_PATH-20);
	strcat(m_sz5UnmergedOutFile,".UPE1");
	strncpy(m_sz3UnmergedOutFile,m_szMergeOutFile,_MAX_PATH-20);
	strcat(m_sz3UnmergedOutFile,".UPE2");
	}
else
	{
	m_sz5UnmergedOutFile[0] = '\0';
	m_sz3UnmergedOutFile[0] = '\0';
	}

m_ProcPhase = ePPRdyOpen;

if((Rslt = ProcOverlapPairs(StartNum)) < eBSFSuccess)
	{
	Reset(false);
	return(Rslt);
	}

Reset(true);
return((int)Rslt);
}

const int cMinReadSeqLen = 20;			// raw read sequences must be at least this length
const int cMaxReadDesrLen = 200;		// allow for read descriptors of no longer than this length

int													// returns number of sequences merged and written to output file
CMergeReadPairs::ProcOverlapPairs(int StartNum)		// initial starting sequence identifier
{
int Rslt;
int Idx;
UINT64 Spacer2 = 0;
UINT8 PE5Seq[cMaxFastQSeqLen+1];
UINT8 PE5Qual[cMaxFastQSeqLen+1];
UINT8 PE3Seq[cMaxFastQSeqLen+1];
UINT8 PE3Qual[cMaxFastQSeqLen+1];
char szPE5DescrBuff[cMaxReadDesrLen+1];
char szPE3DescrBuff[cMaxReadDesrLen+1];
UINT8 MSeq[(2*cMaxFastQSeqLen)+1];
UINT8 szMSeq[(2*cMaxFastQSeqLen)+1];
UINT8 szMQual[(2*cMaxFastQSeqLen)+1];
UINT64 Spacer1 = 0;
int PE5DescrLen;
int PE3DescrLen;
CFasta PE5Fasta;
CFasta PE3Fasta;
int NumPE5Reads;
int NumPE3Reads;
int MinObsOverlap;
int MaxObsOverlap;

int OverlapDistCnts[cMaxFastQSeqLen];
int OverlapDistSubs[cMaxOverlapPercSubs+1];

int MergeIdx;
UINT8 *pMSeq;
UINT8 *pMSeqQual;
UINT8 *pSeq5Qual;
UINT8 *pSeq3Qual;
etSeqBase *pSeq5;
etSeqBase *pSeq3;
etSeqBase Base5;
etSeqBase Base3;
int OL5Idx;
int OL3Idx;
int MaxOverlap;
int ReqOverlap3;
int MaxSubs;
int CurSubs;
int AllowedSubs;

int NumOverlapping;
int PE5SeqLen;
int PE3SeqLen;
int PE5QualLen;
int PE3QualLen;

if((Rslt=(teBSFrsltCodes)PE5Fasta.Open(m_szIn5ReadsFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",m_szIn5ReadsFile,PE5Fasta.ErrText((teBSFrsltCodes)Rslt),PE5Fasta.GetErrMsg());
	Reset(false);
	return(Rslt);
	}

if(PE5Fasta.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"SOLiD colorspace processing is not supported, '%s' is colorspace...",m_szIn5ReadsFile);
	PE5Fasta.Close();
	Reset(false);
	return(eBSFerrOpnFile);
	}

if((Rslt=(teBSFrsltCodes)PE3Fasta.Open(m_szIn3ReadsFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",m_szIn3ReadsFile,PE3Fasta.ErrText((teBSFrsltCodes)Rslt),PE3Fasta.GetErrMsg());
	PE5Fasta.Close();
	Reset(false);
	return(Rslt);
	}

if(PE3Fasta.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"SOLiD colorspace processing is not supported, '%s' is colorspace...",m_szIn3ReadsFile);
	PE5Fasta.Close();
	PE3Fasta.Close();
	Reset(false);
	return(eBSFerrOpnFile);
	}

m_bIsFastq = PE5Fasta.IsFastq();
if(m_bIsFastq != PE3Fasta.IsFastq())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Paired end files are not both either fastq or fasta type\n '%s' is %s, '%s' is %s",
								m_szIn5ReadsFile,m_bIsFastq ? "fastq" : "fasta",  m_szIn3ReadsFile,PE3Fasta.IsFastq() ? "fastq" : "fasta");
	PE5Fasta.Close();
	PE3Fasta.Close();
	Reset(false);
	return(eBSFerrOpnFile);
	}

if(m_OFormat == eOFauto)
	m_OFormat = m_bIsFastq ? eOFfastq : eOFfasta;

if(m_bAppendOut)
	{
#ifdef _WIN32
	m_hOutMerged = open(m_szMergeOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
	m_hOutMerged = open(m_szMergeOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
	if(m_hOutMerged == -1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_szMergeOutFile);
		PE5Fasta.Close();
		PE3Fasta.Close();
		Reset(false);
		return(eBSFerrCreateFile);
		}
	// seek to end of file ready for appending
	lseek(m_hOutMerged,0,SEEK_END);

	if(m_sz5UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
		m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
		if(m_hOut5Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz5UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		// seek to end of file ready for appending
		lseek(m_hOut5Unmerged,0,SEEK_END);
		}

	if(m_sz3UnmergedOutFile[0] != '\0')
		{
#ifdef _WIN32
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
		m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
		if(m_hOut3Unmerged == -1)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output 3' unmerged file '%s'",m_sz3UnmergedOutFile);
			Reset(false);
			return(eBSFerrCreateFile);
			}
		// seek to end of file ready for appending
		lseek(m_hOut3Unmerged,0,SEEK_END);
		}
	}
else
	{
#ifdef _WIN32
	m_hOutMerged = open(m_szMergeOutFile,O_CREATETRUNC );
#else
	if((m_hOutMerged = open(m_szMergeOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutMerged,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_szMergeOutFile,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif
	if(m_hOutMerged == -1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_szMergeOutFile);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	lseek(m_hOutMerged,0,SEEK_SET);

	if(m_PMode == ePMseparate)
		{
		if(m_sz5UnmergedOutFile[0] != '\0')
			{
#ifdef _WIN32
			m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_CREATETRUNC );
#else
			if((m_hOut5Unmerged = open(m_sz5UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
				if(ftruncate(m_hOut5Unmerged,0)!=0)
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_sz5UnmergedOutFile,strerror(errno));
						Reset(false);
						return(eBSFerrCreateFile);
						}
#endif
			if(m_hOut5Unmerged == -1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz5UnmergedOutFile);
				Reset(false);
				return(eBSFerrCreateFile);
				}
			lseek(m_hOut5Unmerged,0,SEEK_SET);
			}

		if(m_sz3UnmergedOutFile[0] != '\0')
			{
#ifdef _WIN32
			m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_CREATETRUNC );
#else
			if((m_hOut3Unmerged = open(m_sz3UnmergedOutFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
				if(ftruncate(m_hOut3Unmerged,0)!=0)
						{
						gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s merged PE microcontigs - %s",m_sz3UnmergedOutFile,strerror(errno));
						Reset(false);
						return(eBSFerrCreateFile);
						}
#endif
			if(m_hOut3Unmerged == -1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output merged PE microcontigs file '%s'",m_sz3UnmergedOutFile);
				Reset(false);
				return(eBSFerrCreateFile);
				}
			lseek(m_hOut3Unmerged,0,SEEK_SET);
			}
		}
	}

if(m_hOutMerged != -1 && m_pszMSeqs == NULL)
	{
	if((m_pszMSeqs = new char [cAllocOutBuffLen])==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d) for buffering file writes",cAllocOutBuffLen);
		Reset(false);
		return(eBSFerrMem);
		}
	m_AllocdMSeqs = cAllocOutBuffLen;
	m_CurMSeqLen = 0;
	}

if(m_PMode == ePMseparate)
	{
	if(m_hOut5Unmerged != -1 && m_pszUnmergedP1Seqs == NULL)
		{
		if((m_pszUnmergedP1Seqs = new char [cAllocOutBuffLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d) for buffering file writes",cAllocOutBuffLen);
			Reset(false);
			return(eBSFerrMem);
			}
		m_AllocdUnmergedP1Seqs = cAllocOutBuffLen;
		m_CurUnmergedP1Seqs = 0;
		}

	if(m_hOut3Unmerged != -1 && m_pszUnmergedP2Seqs == NULL)
		{
		if((m_pszUnmergedP2Seqs = new char [cAllocOutBuffLen])==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to allocate memory (%d) for buffering file writes",cAllocOutBuffLen);
			Reset(false);
			return(eBSFerrMem);
			}
		m_AllocdUnmergedP2Seqs = cAllocOutBuffLen;
		m_CurUnmergedP2Seqs = 0;
		}
	}

m_ProcPhase = ePRdyProc;

memset(OverlapDistCnts,0,sizeof(OverlapDistCnts));
memset(OverlapDistSubs,0,sizeof(OverlapDistSubs));
NumPE5Reads = 0;
NumPE3Reads = 0;
NumOverlapping = 0;
MinObsOverlap = -1;
MaxObsOverlap = -1;
time_t Started = time(0);
while((Rslt = (teBSFrsltCodes)(PE5SeqLen = PE5Fasta.ReadSequence(PE5Seq,cMaxFastQSeqLen,true,false))) > eBSFSuccess)
	{
	if(!(NumPE5Reads % 10000) && NumPE5Reads > 0)
		{
		time_t Now = time(0);
		unsigned long ElapsedSecs = (unsigned long) (Now - Started);
		if(ElapsedSecs >= 60)
			{
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d paired reads of which %d overlapped",NumPE5Reads,NumOverlapping);
			Started = Now;

			if(m_PMode == ePMseparate)
				{
				if(m_CurUnmergedP1Seqs > 0)
					{
					CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,(size_t)m_CurUnmergedP1Seqs);
					m_CurUnmergedP1Seqs = 0;
					}
				if(m_CurUnmergedP2Seqs > 0)
					{
					CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,(size_t)m_CurUnmergedP2Seqs);
					m_CurUnmergedP2Seqs = 0;
					}
				}

			if(m_CurMSeqLen > 0)
				{
				CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
				m_CurMSeqLen = 0;
				}
			}
		}
 
	NumPE5Reads += 1;
	if(PE5SeqLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		PE5DescrLen = PE5Fasta.ReadDescriptor((char *)szPE5DescrBuff,sizeof(szPE5DescrBuff)-1);
		szPE5DescrBuff[sizeof(szPE5DescrBuff)-1] = '\0';
		PE5SeqLen = PE5Fasta.ReadSequence(PE5Seq,cMaxFastQSeqLen);
		if(PE5SeqLen < cMinReadSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence from '%s' after %d reads parsed",m_szIn5ReadsFile,NumPE5Reads);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE5DescrBuff);
			PE5Fasta.Close();
			PE3Fasta.Close();
			return(eBSFerrParse);
			}
		if(m_bIsFastq)
			{
			PE5QualLen = PE5Fasta.ReadQValues((char *)PE5Qual,cMaxFastQSeqLen);
			PE5Qual[PE5QualLen] = '\0';
			if(PE5QualLen != PE5SeqLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing quality from '%s' after %d reads parsed",m_szIn5ReadsFile,NumPE5Reads);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence length: %d Quality length: %d",PE5SeqLen,PE5QualLen);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE5DescrBuff);
				PE5Fasta.Close();
				PE3Fasta.Close();
				return(eBSFerrParse);
				}
			}
		}

	Rslt = (teBSFrsltCodes)(PE3SeqLen = PE3Fasta.ReadSequence(PE3Seq,cMaxFastQSeqLen,true,false));
	if(Rslt <= eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Fewer reads in '%s' than '%s' after %d reads parsed",m_szIn3ReadsFile,m_szIn5ReadsFile,NumPE3Reads);
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE3DescrBuff);
		PE5Fasta.Close();
		PE3Fasta.Close();
		return(eBSFerrParse);
		}	
	NumPE3Reads += 1;
	if(PE3SeqLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		PE3DescrLen = PE3Fasta.ReadDescriptor((char *)szPE3DescrBuff,sizeof(szPE3DescrBuff)-1);
		szPE3DescrBuff[sizeof(szPE3DescrBuff)-1] = '\0'; 
		PE3SeqLen = PE3Fasta.ReadSequence(PE3Seq,cMaxFastQSeqLen);
		if(PE3SeqLen < cMinReadSeqLen)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence from '%s' after %d reads parsed",m_szIn3ReadsFile,NumPE3Reads);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE3DescrBuff);
			PE5Fasta.Close();
			PE3Fasta.Close();
			return(eBSFerrParse);
			}
		if(m_bIsFastq)
			{
			PE3QualLen = PE3Fasta.ReadQValues((char *)PE3Qual,cMaxFastQSeqLen);
			PE3Qual[PE3QualLen] = '\0';
			if(PE3QualLen != PE3SeqLen)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing quality from '%s' after %d reads parsed",m_szIn3ReadsFile,NumPE3Reads);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequence length: %d Quality length: %d",PE3SeqLen,PE3QualLen);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szPE5DescrBuff);
				PE5Fasta.Close();
				PE3Fasta.Close();
				return(eBSFerrParse);
				}
			}
		}

	// PE3 needs to be revcpl'd
	CSeqTrans::ReverseComplement(PE3SeqLen,PE3Seq);
	if(m_bIsFastq)
		CSeqTrans::ReverseSeq(PE3QualLen,PE3Qual);

	// now try for maximal overlap of at least m_MinOverlap allowing (if user specified) for sequencer base call errors 
	MaxSubs = m_MaxOverlapPropSubs; 
	MaxOverlap = 0;
	for(OL5Idx = 0; OL5Idx <= (PE5SeqLen - m_MinOverlap); OL5Idx++)
		{
		pSeq5 = &PE5Seq[OL5Idx];
		pSeq3 = PE3Seq;
		CurSubs = 0;
		ReqOverlap3 = min(PE5SeqLen - OL5Idx,PE3SeqLen);
		if(MaxSubs != 0)
			{
			if(ReqOverlap3 >= 20)
				AllowedSubs = min(MaxSubs,1 + ((ReqOverlap3 * m_MaxOverlapPropSubs) / 100));
			else
				{
				if(ReqOverlap3 >= 10)
					AllowedSubs = 2;
				else
					{
					if(ReqOverlap3 >= 5)
						AllowedSubs = 1;
					else
						AllowedSubs = 0;
					}
				AllowedSubs = min(MaxSubs,AllowedSubs);
				}
			}
		else
			AllowedSubs = 0;
		for(OL3Idx=0; OL3Idx < ReqOverlap3 && CurSubs <= AllowedSubs; OL3Idx++,pSeq5++,pSeq3++)
			{
			Base5 = *pSeq5 & 0x07;
			Base3 = *pSeq3 & 0x07;
			if(Base5 > eBaseT || Base3 > eBaseT || Base5 != Base3)
				{
				CurSubs += 1;
				if(CurSubs > AllowedSubs)
					break;
				}
			}

		if(OL3Idx == ReqOverlap3 && CurSubs <= MaxSubs)
			{
			if(CurSubs == 0 || CurSubs < MaxSubs)
				MaxOverlap = OL3Idx;
			MaxSubs = CurSubs;
			if(MaxSubs == 0)
				break;
			}
		}

	if(MaxOverlap < 1 && m_PMode != ePMdefault)
		{
		CSeqTrans::ReverseComplement(PE3SeqLen,PE3Seq);
		CSeqTrans::MapSeq2UCAscii(PE3Seq,PE3SeqLen,(char *)szMSeq);
		szMSeq[PE3SeqLen] = '\0';
		if(m_OFormat == eOFfastq)
			{
			if(m_bIsFastq)
				CSeqTrans::ReverseSeq(PE3QualLen,PE3Qual);
			else
				memset(PE3Qual,cPhredOrphanScore,PE3SeqLen);
			PE3Qual[PE3SeqLen] = '\0';
			}

		if(m_PMode == ePMseparate)
			{
			if(m_OFormat == eOFfasta)
				m_CurUnmergedP2Seqs += sprintf((char *)&m_pszUnmergedP2Seqs[m_CurUnmergedP2Seqs],">%s\n%s\n",szPE3DescrBuff,(char *)szMSeq);
			else
				m_CurUnmergedP2Seqs += sprintf((char *)&m_pszUnmergedP2Seqs[m_CurUnmergedP2Seqs],"@%s\n%s\n+\n%s\n",szPE3DescrBuff,(char *)szMSeq,(char *)PE3Qual);
				
			if(m_CurUnmergedP2Seqs > (cAllocOutBuffLen - (8 * cMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,(size_t)m_CurUnmergedP2Seqs);
				m_CurUnmergedP2Seqs = 0;
				}
			}
		else
			{
			if(m_OFormat == eOFfasta)
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],">%s\n%s\n",szPE3DescrBuff,(char *)szMSeq);
			else
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],"@%s\n%s\n+\n%s\n",szPE3DescrBuff,(char *)szMSeq,(char *)PE3Qual);

			if(m_CurMSeqLen > (cAllocOutBuffLen - (8 * cMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
				m_CurMSeqLen = 0;
				}	
			}

		CSeqTrans::MapSeq2UCAscii(PE5Seq,PE5SeqLen,(char *)szMSeq);
		szMSeq[PE5SeqLen] = '\0';
		if(m_OFormat == eOFfastq)
			{
			if(!m_bIsFastq)
				memset(PE5Qual,cPhredOrphanScore,PE5SeqLen);
			PE5Qual[PE5SeqLen] = '\0';
			}

			
		if(m_PMode == ePMseparate)
			{
			if(m_OFormat == eOFfasta)
				m_CurUnmergedP1Seqs += sprintf((char *)&m_pszUnmergedP1Seqs[m_CurUnmergedP1Seqs],">%s\n%s\n",szPE5DescrBuff,(char *)szMSeq);
			else
				m_CurUnmergedP1Seqs += sprintf((char *)&m_pszUnmergedP1Seqs[m_CurUnmergedP1Seqs],"@%s\n%s\n+\n%s\n",szPE5DescrBuff,(char *)szMSeq,(char *)PE5Qual);

			if(m_CurUnmergedP1Seqs > (cAllocOutBuffLen - (8 * cMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,(size_t)m_CurUnmergedP1Seqs);
				m_CurUnmergedP1Seqs = 0;
				}
			}
		else
			{
			if(m_OFormat == eOFfasta)
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],">%s\n%s\n",szPE5DescrBuff,(char *)szMSeq);
			else
				m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],"@%s\n%s\n+\n%s\n",szPE5DescrBuff,(char *)szMSeq,(char *)PE5Qual);

			if(m_CurMSeqLen > (cAllocOutBuffLen - (8 * cMaxFastQSeqLen)))
				{
				CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
				m_CurMSeqLen = 0;
				}	
			}

		continue;
		}

	if(MaxOverlap > 0)
		{
		if(MinObsOverlap == -1 || MaxOverlap < MinObsOverlap)
			MinObsOverlap = MaxOverlap;
		if(MaxObsOverlap == -1 || MaxOverlap > MaxObsOverlap)
			MaxObsOverlap = MaxOverlap;
		OverlapDistSubs[MaxSubs] += 1;
		OverlapDistCnts[MaxOverlap] += 1;
		NumOverlapping += 1;

		// output to file...
		pMSeq = MSeq;
		pSeq5 = PE5Seq;
		pMSeqQual = szMQual;
		pSeq5Qual = PE5Qual;
		for(MergeIdx = 0; MergeIdx < OL5Idx; MergeIdx++,pSeq5++,pMSeq++,pMSeqQual++,pSeq5Qual++)
			{
			*pMSeq = *pSeq5;
			if(m_OFormat == eOFfastq)
				{
				if(m_bIsFastq)
					*pMSeqQual = *pSeq5Qual;
				else
					*pMSeqQual = cPhredHiScore;
				}
			}
		
		pSeq3 = PE3Seq;
		pSeq3Qual = PE3Qual;

		if(MergeIdx < PE5SeqLen)
			for(; MergeIdx < PE5SeqLen; MergeIdx++,pSeq5++,pMSeq++,pSeq3++,pMSeqQual++,pSeq5Qual++,pSeq3Qual++)
				{
				if(*pSeq5 == *pSeq3)
					{
					*pMSeq = *pSeq3;
					if(m_OFormat == eOFfastq)
						{
						if(m_bIsFastq)
							{
							if(*pSeq3Qual >= *pSeq5Qual)
								*pMSeqQual = *pSeq3Qual;
							else
								*pMSeqQual = *pSeq5Qual;
							}
						else
							*pMSeqQual = cPhredHiScore;
						}
					}
				else	// base difference: wheich one to choose... 
					{
					if(m_bIsFastq)	// choose base with highest Phred score
						{
						if(*pSeq3Qual >= *pSeq5Qual)
							{
							*pMSeq = *pSeq3;
							if(m_OFormat == eOFfastq)
								*pMSeqQual = *pSeq3Qual;
							}
						else
							{
							*pMSeq = *pSeq5;
							if(m_OFormat == eOFfastq)
								*pMSeqQual = *pSeq5Qual;
							}
						}
					else
						{
						*pMSeq = *pSeq3;
						if(m_OFormat == eOFfastq)
							*pMSeqQual = cPhredLowScore;
						}
					}
				}

			MergeIdx = PE3SeqLen - OL5Idx;
			if(MergeIdx < PE3SeqLen)
				for( ;MergeIdx < PE3SeqLen;MergeIdx++,pMSeq++,pSeq3++,pMSeqQual++,pSeq3Qual++)
					{
					*pMSeq = *pSeq3;
					if(m_OFormat == eOFfastq)
						{
						if(m_bIsFastq)
							*pMSeqQual = *pSeq3Qual;
						else
							*pMSeqQual = cPhredHiScore;
						}
					}

		CSeqTrans::MapSeq2UCAscii(MSeq,OL5Idx+PE3SeqLen,(char *)szMSeq);
		szMSeq[OL5Idx+PE3SeqLen] = '\0';
		if(m_OFormat == eOFfasta)
			m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],">MSeq%d %s\n%s\n",StartNum+NumOverlapping,szPE5DescrBuff,(char *)szMSeq);
		else
			{
			szMQual[OL5Idx+PE3SeqLen] = '\0';
			m_CurMSeqLen += sprintf((char *)&m_pszMSeqs[m_CurMSeqLen],"@MSeq%d %s\n%s\n+\n%s\n",StartNum+NumOverlapping,szPE5DescrBuff,(char *)szMSeq,szMQual);
			}
		if(m_CurMSeqLen > (cAllocOutBuffLen - (8 * cMaxFastQSeqLen)))
			{
			CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
			m_CurMSeqLen = 0;
			}
		}
	}

if(m_CurMSeqLen > 0)
	{
	CUtility::SafeWrite(m_hOutMerged,m_pszMSeqs,(size_t)m_CurMSeqLen);
	m_CurMSeqLen = 0;
	}

if(m_CurUnmergedP1Seqs > 0)
	{
	CUtility::SafeWrite(m_hOut5Unmerged,m_pszUnmergedP1Seqs,(size_t)m_CurUnmergedP1Seqs);
	m_CurUnmergedP1Seqs = 0;
	}

if(m_CurUnmergedP2Seqs > 0)
	{
	CUtility::SafeWrite(m_hOut3Unmerged,m_pszUnmergedP2Seqs,(size_t)m_CurUnmergedP2Seqs);
	m_CurUnmergedP2Seqs = 0;
	}


gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processed %d paired reads of which %d overlapped",NumPE5Reads,NumOverlapping);
if(NumOverlapping)
	{
	for(Idx = MinObsOverlap; Idx <= MaxObsOverlap; Idx++)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"OverlapLen: %d Counts: %d",Idx,OverlapDistCnts[Idx]);
	for(Idx = 0; Idx <= cMaxOverlapPercSubs; Idx++)
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Substitutions: %d Counts: %d",Idx,OverlapDistSubs[Idx]);
	}

PE5Fasta.Close();
PE3Fasta.Close();

return(NumOverlapping);
}
