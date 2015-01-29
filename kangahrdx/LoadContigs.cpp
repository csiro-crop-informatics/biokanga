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

#include "HomozyReduce.h"

static UINT8 *m_xpConcatSeqs;

teBSFrsltCodes
CHomozyReduce::Disk2Hdr(tsBSFRdsHdr *pRdsHeader,	// load header into this
						char *pszRdsFile)			// loading is from this file
{
if(_lseeki64(m_hInFile,0,SEEK_SET)!=0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Seek failed to offset 0 - %s",pszRdsFile,strerror(errno));
	Reset(false);			// closes opened file..
	return(eBSFerrFileAccess);
	}

// read in header..
if(sizeof(tsBSFRdsHdr) != read(m_hInFile,pRdsHeader,sizeof(tsBSFRdsHdr)))
	{
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrNotBioseq);
	}

// header read, validate it as being a reads file header
if(tolower(pRdsHeader->Magic[0]) != 'b' ||
	tolower(pRdsHeader->Magic[1]) != 'i' ||
	tolower(pRdsHeader->Magic[2]) != 'o' ||
	tolower(pRdsHeader->Magic[3]) != 'r')
	{
	close(m_hInFile);
	m_hInFile = -1;
	return(eBSFerrNotBioseq);
	}

	// can we handle this version?
if(pRdsHeader->Version < cBSFRRRdsVersionBack || pRdsHeader->Version > cBSFRRRdsVersion)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"%s opened as a pre-processed reads file - expected between version %d and %d, file version is %d",
					pszRdsFile,cBSFRRRdsVersionBack,cBSFRRRdsVersion,pRdsHeader->Version);
	Reset(false);			// closes opened file..
	return(eBSFerrFileVer);
	}

return(eBSFSuccess);
}


// LoadRawContigs
// Load reads from fasta, fastq or csfasta formated raw reads file
teBSFrsltCodes
CHomozyReduce::LoadRawContigs(int MaxNs,		// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					int FileID,				// uniquely identifies source file
					char *pszFile)			// process from this file
{
static int FileNamesOfs = 0;
teBSFrsltCodes Rslt;
int NumDescrContigs;
int NumUnderlen;
int NumAcceptedContigs;
bool bIsFastq;
int ReadLen;

int DescrLen;
UINT8 szDescrBuff[1024];

int NumInvalValues = 0;
int NumUnsupportedBases = 0;
int NumUnderlength = 0;
int NumExcessNs = 0;

etSeqBase *pSeq;
UINT8 *pRawContigsBuff;
if((pRawContigsBuff = new UINT8 [cAllocRawSeqLen]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for raw reads buffering...");
	Reset(false);
	return(eBSFerrMem);
	}

CFasta Fasta;
if((Rslt=(teBSFrsltCodes)Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	delete pRawContigsBuff;
	Reset(false);
	return(Rslt);
	}

if(Fasta.IsSOLiD())
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Sequences from '%s' are in SOLiD colorspace...",pszFile);
	delete pRawContigsBuff;
	Reset(false);
	return(eBSFerrOpnFile);
	}

bIsFastq = Fasta.IsFastq();

// get file size and alloc for expected total sequences up front to try and reduce the number of reallocs which are required
// may end up allocating more memory than actually needed but on Windows HPC seems to result in a significant throughput improvement
UINT64 ReqAllocSize = Fasta.InitialFileSize();
if(bIsFastq)
	ReqAllocSize /= 2;			// quality scores effectively double the file size relative to the actual sequences

AcquireSerialise();
if(m_pSeqs2Assemb == NULL)
	{
	m_AllocMemSeqs2Assemb = (size_t)ReqAllocSize;
#ifdef _WIN32
	m_pSeqs2Assemb = (UINT8 *) malloc((size_t)m_AllocMemSeqs2Assemb);
	if(m_pSeqs2Assemb == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Concatenated sequences memory allocation of %lld bytes - %s",(INT64)m_AllocMemSeqs2Assemb,strerror(errno));
		m_AllocMemConcat = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqs2Assemb = (UINT8 *)mmap(NULL,m_AllocMemSeqs2Assemb, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqs2Assemb == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Concatenated sequences memory of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemSeqs2Assemb,strerror(errno));
		m_pSeqs2Assemb = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_Seqs2AssembLen = 0;
	m_NumSeqs2Assemb = 0;
	}
else
	{
	UINT8 *pDstSeq;
	size_t memreq;
	if((m_Seqs2AssembLen + ReqAllocSize + 100) >= m_AllocMemSeqs2Assemb)		// 100 as a small safety margin!
		{
		memreq = (size_t)(m_Seqs2AssembLen + ReqAllocSize);
#ifdef _WIN32
		pDstSeq = (UINT8 *) realloc(m_pSeqs2Assemb,memreq);
#else
		pDstSeq = (UINT8 *)mremap(m_pSeqs2Assemb,m_AllocMemSeqs2Assemb,memreq,MREMAP_MAYMOVE);
		if(pDstSeq == MAP_FAILED)
			pDstSeq = NULL;
#endif
		if(pDstSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocMemSeqs2Assemb = memreq;
		m_pSeqs2Assemb = pDstSeq;
		}
	}
ReleaseSerialise();

if(m_AllocMemSeqs2Assemb > m_CurMaxMemWorkSetBytes)
	{
	if(!SetMaxMemWorkSetSize(m_AllocMemSeqs2Assemb * 2))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Length of loaded concatenated sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
		return(eBSFerrMaxDirEls);
		}
	}

m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles].SrcFileID = FileID;
strncpy((char *)m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles].SrcFileName,pszFile,sizeof(m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles].SrcFileName)-1);

NumUnsupportedBases = 0;
NumDescrContigs = 0;
NumUnderlen = 0;
while((Rslt = (teBSFrsltCodes)(ReadLen = Fasta.ReadSequence(pRawContigsBuff,cAllocRawSeqLen,true,false))) > eBSFSuccess)
	{
	NumDescrContigs += 1;
	if(ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		DescrLen = Fasta.ReadDescriptor((char *)szDescrBuff,sizeof(szDescrBuff)-1);
		szDescrBuff[sizeof(szDescrBuff)-1] = '\0';
		ReadLen = Fasta.ReadSequence(pRawContigsBuff,cAllocRawSeqLen);
		if(ReadLen < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed",NumDescrContigs);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szDescrBuff);
			delete pRawContigsBuff;
			Fasta.Close();
			return(eBSFerrParse);
			}

		if(ReadLen < (Trim5 + Trim3 + MinSeqLen))
			{
			NumUnderlen += 1;
			continue;
			}

		pSeq = pRawContigsBuff;

		// trim 5' and 3' as requested
		if(Trim5 > 0)
			{
			pSeq += Trim5;
			ReadLen -= Trim5;
			}

		ReadLen -= Trim3;

		// check for excessive number of Ns
		int Idx;
		int NumNs = 0;		// number of indeterminate bases in last 100bp window
		etSeqBase *pBase = pSeq;
		for(Idx = 0; Idx < ReadLen; Idx++,pBase++)
			{
			if(Idx >= 100)
				{
				if((pBase[-100] & 0x07) == eBaseN)
					NumNs -= 1;
				}
			if((*pBase & 0x07) == eBaseN)
				NumNs += 1;
			if(NumNs > MaxNs)
				break;
			}

		if(NumNs > MaxNs ||
				(ReadLen <= 100 && NumNs > ((MaxNs * ReadLen)/100)))
			{
			NumExcessNs += 1;
			continue;
			}

		if((Rslt=AddSeq(ReadLen,pSeq)) !=eBSFSuccess)
			break;

		if(m_Seqs2AssembLen > cMaxConcatSeqLen)			// hit limit?
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Length of loaded concatenated sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
			Rslt = eBSFerrMaxDirEls;
			break;
			}
		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszFile,
				Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
		Fasta.Close();
		return(eBSFerrParse);
		}
	}
if(pRawContigsBuff != NULL)
	delete pRawContigsBuff;

if(Rslt != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszFile);
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	Fasta.Close();
	return(Rslt);
	}
Fasta.Close();

NumAcceptedContigs = NumDescrContigs - (NumUnderlen + NumExcessNs);
m_TotSeqsParsed += NumDescrContigs;
m_TotSeqsUnderLen += NumUnderlen;
m_TotSeqsExcessNs += NumExcessNs;
m_TotSeqs2Assemb += NumAcceptedContigs;
m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles++].NumContigs = NumAcceptedContigs;
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContigs: Completed parsing %d sequences from '%s'",NumDescrContigs,pszFile);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContigs: %d were under length %d, %d excessive indeterminate (%d Ns), sequences retained: %1.9d",NumUnderlen,MinSeqLen,NumExcessNs,MaxNs, NumAcceptedContigs);
m_NumRawFiles += 1;
return((teBSFrsltCodes)NumAcceptedContigs);
}

// LoadContigs
// Load pre-processed reads (as generated by kangar) from '.rds' structured file
// If unable to load as .rds then will try and load with LoadRawContigs()
teBSFrsltCodes
CHomozyReduce::LoadContigs(int MaxNs,				// filter out input sequences having higher than this number of indeterminate bases per 100bp (default is 1, range 0..10)
					int Trim5,				// trim this number of 5' bases from input sequences (default is 0, range 0..20)
					int Trim3,				// trim this number of 3' bases from input sequences (default is 0, range 0..20)
					int MinSeqLen,		    // filter out input sequences (after any trimming) which are less than this length (default is 50bp, range 30..10000)
					 int FileID,			// uniquely identifies source file
					 char *pszFile)			// file containing contigs
{
teBSFrsltCodes Rslt;
int NumContigs;
int NumUnderlen;
int NumExcessNs;
int RdLen;
tsRawReadV5 *pReadV5;						// current pre-processed read being processed
tsRawReadV6 *pReadV6;						// current pre-processed read being processed
int SizeOfRawRead;
int CurReadLen;
int CurDescrLen;

UINT8 ReadBuff[0x0ffff];				// used to buffer pre-processed reads

etSeqBase *pSeq;
int ReadLen;
int NumAcceptedContigs;

int BuffLen;
int BuffOfs;
UINT8 *pSrcSeq;							// used when copying read sequence (src) into concatenated sequence (dst)
tsBSFRdsHdr RdsHdr;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading sequences from: %s",pszFile);

#ifdef _WIN32
m_hInFile = open(pszFile, O_READSEQ ); // file access is normally sequential..
#else
m_hInFile = open64(pszFile, O_READSEQ );
#endif

if(m_hInFile == -1)					// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszFile,strerror(errno));
	Reset(false);
	return(eBSFerrOpnFile);
	}

// expecting a pre-processed .rds file as input, header processing will confirm!
// if not a bioseq .rds then will try processing as fasta
if((Rslt=Disk2Hdr(&RdsHdr,pszFile))!=eBSFSuccess)
	{
	if(Rslt == eBSFerrNotBioseq)
		return(LoadRawContigs(MaxNs,Trim5,Trim3,MinSeqLen,FileID,pszFile));		// try as raw fasta
	Reset(true);
	return((teBSFrsltCodes)Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contigs file '%s' generator was version: %d",pszFile,RdsHdr.Version);
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contains %d sequences from %d source sequences with duplicate sequences %s, PMode was %d, %d bases 5' trimmed, %d bases 3' trimmed",
		RdsHdr.NumRds,RdsHdr.OrigNumReads,RdsHdr.FlagsK ? "retained":"removed", RdsHdr.PMode,RdsHdr.Trim5,RdsHdr.Trim3);
gDiagnostics.DiagOutMsgOnly(eDLInfo,"Sequences were originally processed from %d files",RdsHdr.NumFiles);
char *pszSrcFile = (char *)RdsHdr.FileNames;
for(BuffOfs=0;BuffOfs<RdsHdr.NumFiles;BuffOfs++)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source file: '%s'",pszSrcFile);
	pszSrcFile += strlen(pszSrcFile) + 1;
	}

// ensure there is at least one read...
if(RdsHdr.NumRds == 0 || RdsHdr.TotReadsLen == 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Nothing to do, '%s' contains no reads...",pszFile);
	Reset(true);
	return(eBSFerrOpnFile);
	}

if(RdsHdr.FlagsCS)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Basespace assembly requested but sequences from '%s' are in SOLiD colorspace...",pszFile);
	Reset(false);
	return(eBSFerrOpnFile);
	}

m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles].SrcFileID = FileID;
strncpy((char *)m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles].SrcFileName,pszFile,sizeof(m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles].SrcFileName)-1);

BuffLen = 0;
BuffOfs = 0;
NumContigs = 0;
NumUnderlen = 0;
NumExcessNs = 0;

// tsRawRead struct size is version dependent
SizeOfRawRead = RdsHdr.Version < 6 ? sizeof(tsRawReadV5) : sizeof(tsRawReadV6);

// iterate each read sequence starting from the first
lseek(m_hInFile,(long)RdsHdr.RdsOfs,SEEK_SET);
while((RdLen = read(m_hInFile,&ReadBuff[BuffLen],sizeof(ReadBuff) - BuffLen)) > 0)
	{
	BuffLen += RdLen;
	BuffOfs = 0;
	while((BuffLen - BuffOfs) >=  SizeOfRawRead)
		{
		if(RdsHdr.Version < 6)
			{
			pReadV5 = (tsRawReadV5 *)&ReadBuff[BuffOfs];
			pSrcSeq = &pReadV5->Read[pReadV5->DescrLen+1];
			CurDescrLen = pReadV5->DescrLen;
			CurReadLen = pReadV5->ReadLen;
			}
		else
			{
			pReadV6 = (tsRawReadV6 *)&ReadBuff[BuffOfs];
			pSrcSeq = &pReadV6->Read[pReadV6->DescrLen+1];
			CurDescrLen = pReadV6->DescrLen;
			CurReadLen = pReadV6->ReadLen;
			}


		if((int)(CurDescrLen + CurReadLen + SizeOfRawRead) > (BuffLen - BuffOfs))
			break;
		BuffOfs += CurDescrLen + CurReadLen + SizeOfRawRead;
		// pRead now pts at a read

		NumContigs += 1;
		pSeq = pSrcSeq;
		ReadLen = CurReadLen;

		if(ReadLen < (Trim5 + Trim3 + MinSeqLen))
			{
			NumUnderlen += 1;
			continue;
			}

		// trim 5' and 3' as requested
		if(Trim5 > 0)
			{
			pSeq += Trim5;
			ReadLen -= Trim5;
			}
		ReadLen -= Trim3;

		// check for excessive number of Ns
		int Idx;
		int NumNs = 0;		// number of indeterminate bases in last 100bp window
		etSeqBase *pBase = pSeq;
		for(Idx = 0; Idx < ReadLen; Idx++,pBase++)
			{
			if(Idx >= 100)
				{
				if((pBase[-100] & 0x07) == eBaseN)
					NumNs -= 1;
				}
			if((*pBase & 0x07) == eBaseN)
				NumNs += 1;
			if(NumNs > MaxNs)
				break;
			}

		if(NumNs > MaxNs ||
				(ReadLen <= 100 && NumNs > ((MaxNs * ReadLen)/100)))
			{
			NumExcessNs += 1;
			continue;
			}

		if(m_Seqs2AssembLen * 2 > m_CurMaxMemWorkSetBytes)
			{
			if(!SetMaxMemWorkSetSize(m_Seqs2AssembLen * 5))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Length of loaded concatenated sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
				Rslt = eBSFerrMaxDirEls;
				break;
				}
			}
		if((Rslt=AddSeq(ReadLen,pSeq))!=eBSFSuccess)
				break;

		}
	if(m_RdsSfxHdr.ConcatSeqLen > cMaxConcatSeqLen)			// hit limit?
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Loaded length of concatenated sequences exceeds limit of %lld bytes",cMaxConcatSeqLen);
		Rslt = eBSFerrMaxDirEls;
		break;
		}

	BuffLen -= BuffOfs;
	if(BuffLen)
		memcpy(ReadBuff,&ReadBuff[BuffOfs],BuffLen);
	}
close(m_hInFile);		// reads all loaded from the .rds file so can now close
m_hInFile = -1;

NumAcceptedContigs = NumContigs - (NumUnderlen + NumExcessNs);
m_TotSeqsParsed += NumContigs;
m_TotSeqsUnderLen += NumUnderlen;
m_TotSeqsExcessNs += NumExcessNs;
m_TotSeqs2Assemb += NumAcceptedContigs;
m_RdsSfxHdr.RdsSrcFiles[m_RdsSfxHdr.NumSrcFiles++].NumContigs = NumAcceptedContigs;
if(Rslt == eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContigs: Completed parsing %d sequences from '%s'",NumContigs,pszFile);
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"LoadContigs: %d were under length %d, %d excessive indeterminate (%d Ns), sequences retained: %1.9d",NumUnderlen,MinSeqLen,NumExcessNs,MaxNs,NumAcceptedContigs);
	}
m_NumRawFiles += 1;
return(Rslt == eBSFSuccess ? (teBSFrsltCodes)NumAcceptedContigs : Rslt);
}

// GenSfxdSeqs
// Assumes that on entry m_pSeqs2Assemb holds m_NumSeqs2Assemb and is of m_Seqs2AssembLen
// Writes concatenated sequences into m_pConcatSeqs
teBSFrsltCodes								// returns number of concatenated sequences or if < 0 then error code
CHomozyReduce::GenSfxdSeqs(void)			// generates concatenated sequences into m_pConcatSeqs with m_pSfxdSeqs holding offsets into m_pConcatSeqs ready for sorting
{
int Idx;
size_t MemReq;
tsSfxdSeq *pRead;				// current read
UINT8 *pDstSeq;
UINT64 XFormID;
UINT32 SeqIdx;
UINT8 *pSeq;
UINT8 *pNxtSeq;
UINT32 SeqLen;

m_NumSfxdSeqs = 0;
m_RdsSfxHdr.SumReadLens = 0;
m_RdsSfxHdr.ConcatSeqLen = 0;

if(m_pSeqs2Assemb == NULL || m_Seqs2AssembLen == 0 || m_NumSeqs2Assemb == 0)
	return((teBSFrsltCodes)0);

if(m_pSfxdSeqs == NULL)					// if first time then memory needs to be allocated
	{
	m_AllocdNumSfxdSeqs = m_NumSeqs2Assemb + 10;		// small safety margin
	m_AllocdMemSfxdSeqs = ((UINT64)m_AllocdNumSfxdSeqs * sizeof(tsSfxdSeq));
#ifdef _WIN32
	m_pSfxdSeqs = (tsSfxdSeq *) malloc((size_t)m_AllocdMemSfxdSeqs);
	if(m_pSfxdSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSfxdSeqs: Initial reads memory allocation of %lld bytes - %s",(INT64)m_AllocdMemSfxdSeqs,strerror(errno));
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSfxdSeqs = (tsSfxdSeq *)mmap(NULL,m_AllocdMemSfxdSeqs, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSfxdSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Initial reads memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_AllocdMemSfxdSeqs,strerror(errno));
		m_pSfxdSeqs = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	}
else
	{
	if(m_NumSeqs2Assemb >=  m_AllocdNumSfxdSeqs || ((m_NumSfxdSeqs * 11) < (m_AllocdNumSfxdSeqs * 10)))
		{
		MemReq = (size_t)((UINT64)(m_NumSeqs2Assemb + 10) * sizeof(tsSfxdSeq));
	#ifdef _WIN32
		pRead = (tsSfxdSeq *) realloc(m_pSfxdSeqs,(size_t)MemReq);
	#else
		pRead = (tsSfxdSeq *)mremap(m_pSfxdSeqs,m_AllocdMemSfxdSeqs,MemReq,MREMAP_MAYMOVE);
		if(pRead == MAP_FAILED)
			pRead = NULL;
	#endif
		if(pRead == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSfxdSeqs: Memory re-allocation to %d bytes - %s",MemReq,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocdMemSfxdSeqs = MemReq;
		m_AllocdNumSfxdSeqs = m_NumSeqs2Assemb;
		m_pSfxdSeqs = pRead;
		}
	}

MemReq = (size_t)((UINT64)m_NumSeqs2Assemb * 6) + m_Seqs2AssembLen;
if(m_pConcatSeqs == NULL)
	{
	m_AllocMemConcat = MemReq;
#ifdef _WIN32
	m_pConcatSeqs = (UINT8 *) malloc((size_t)m_AllocMemConcat);
	if(m_pConcatSeqs == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSfxdSeqs: Concatenated sequences memory allocation of %lld bytes - %s",(INT64)m_AllocMemConcat,strerror(errno));
		m_AllocMemConcat = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pConcatSeqs = (UINT8 *)mmap(NULL,m_AllocMemConcat, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pConcatSeqs == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSfxdSeqs: Concatenated sequences memory of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemConcat,strerror(errno));
		m_pConcatSeqs = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	}
else
	{
	// need to account for sequence separator (6 bytes) plus 4 trailing cCSeqEOS plus a few spare
	if(MemReq >= m_AllocMemConcat || (MemReq * 11 < m_AllocMemConcat * 10))
		{
#ifdef _WIN32
		pDstSeq = (UINT8 *) realloc(m_pConcatSeqs,MemReq);
#else
		pDstSeq = (UINT8 *)mremap(m_pConcatSeqs,m_AllocMemConcat,MemReq,MREMAP_MAYMOVE);
		if(pDstSeq == MAP_FAILED)
			pDstSeq = NULL;
#endif
		if(pDstSeq == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenSfxdSeqs: Memory re-allocation to %d bytes - %s",MemReq,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocMemConcat = MemReq;
		m_pConcatSeqs = pDstSeq;
		}
	}

SetMaxMemWorkSetSize((size_t)(m_AllocMemConcat + m_AllocMemSeqs2Assemb + m_AllocMemSfx));

pDstSeq = m_pConcatSeqs;

pNxtSeq = m_pSeqs2Assemb;
for(SeqIdx = 0; SeqIdx < m_NumSeqs2Assemb; SeqIdx++)
	{
	pSeq = pNxtSeq;
	SeqLen = *(UINT32 *)pSeq;
	pSeq += sizeof(UINT32);
	pNxtSeq = pSeq + SeqLen;

	pRead = &m_pSfxdSeqs[m_NumSfxdSeqs++];
	pRead->Flags = 0;
	pRead->SeqID = m_NumSfxdSeqs;

	if(m_RdsSfxHdr.ConcatSeqLen > 1)		// overwrite previously appended cCSeqEOS
		m_RdsSfxHdr.ConcatSeqLen -= 1;
	pDstSeq = &m_pConcatSeqs[m_RdsSfxHdr.ConcatSeqLen];
	if(m_RdsSfxHdr.ConcatSeqLen == 0)
		*pDstSeq++ = cCSeqBOS;
	else
		*pDstSeq++ = cCSeqSep;
	m_RdsSfxHdr.ConcatSeqLen += 1;
	pRead->ConcatSeqOfs = m_RdsSfxHdr.ConcatSeqLen;
	XFormID = IDtoXForm(pRead->SeqID);
	*(UINT32 *)pDstSeq = (UINT32)XFormID;
	pDstSeq += sizeof(UINT32);
	*pDstSeq++ = (UINT8)(XFormID >> 32);
	m_RdsSfxHdr.ConcatSeqLen += 5;

	for(Idx = 0; Idx < (int)SeqLen; Idx++)
		*pDstSeq++ = *pSeq++ & 0x07;
	pRead->ReadLen = SeqLen;

	m_RdsSfxHdr.ConcatSeqLen += pRead->ReadLen;
	*pDstSeq++ = cCSeqEOS;		// will be overwritten by next AddEntry()
	m_RdsSfxHdr.ConcatSeqLen += 1;
	m_RdsSfxHdr.SumReadLens += pRead->ReadLen;

	if(m_RdsSfxHdr.ConcatSeqLen > cMaxConcatSeqLen)
		{
		gDiagnostics.DiagOut(eDLWarn,gszProcName,"GenSfxdSeqs: Reached concatenated sequence size limit of %lld index elements",cMaxConcatSeqLen);
		return(eBSFerrMaxEntries);
		}
	}
m_MeanReadLen = (int)(m_RdsSfxHdr.SumReadLens / m_NumSeqs2Assemb);

m_Seqs2AssembLen = 0;
m_NumSeqs2Assemb = 0;
m_NumSeqsProc = 0;

return((teBSFrsltCodes)m_NumSfxdSeqs);
}

// AddSeq
// Add sequences which are to be subsequently suffix indexed and assembled
teBSFrsltCodes
CHomozyReduce::AddSeq(int SeqLen,				// sequence length
						UINT8 *pSeq)			// ptr to sequence
{
int Idx;
size_t memreq;
UINT8 *pDstSeq;
if(SeqLen < 1 || pSeq == NULL)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"AddSeq: SeqLen: %d pSeq is %s",SeqLen,pSeq == NULL ? "NULL" : "none-null");
	return(eBSFerrParams);
	}

AcquireSerialise();
if(m_pSeqs2Assemb == NULL)
	{
	m_AllocMemSeqs2Assemb = (size_t)cReallocConcatSeqs;
#ifdef _WIN32
	m_pSeqs2Assemb = (UINT8 *) malloc((size_t)m_AllocMemSeqs2Assemb);
	if(m_pSeqs2Assemb == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Concatenated sequences memory allocation of %lld bytes - %s",(INT64)m_AllocMemSeqs2Assemb,strerror(errno));
		m_AllocMemConcat = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSeqs2Assemb = (UINT8 *)mmap(NULL,m_AllocMemSeqs2Assemb, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSeqs2Assemb == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadContigs: Concatenated sequences memory of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemSeqs2Assemb,strerror(errno));
		m_pSeqs2Assemb = NULL;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	m_Seqs2AssembLen = 0;
	m_NumSeqs2Assemb = 0;
	}

if((m_Seqs2AssembLen + SeqLen + 100) >= m_AllocMemSeqs2Assemb)		// 100 as a small safety margin!
	{
	memreq = (size_t)(m_Seqs2AssembLen + cReallocConcatSeqs);
#ifdef _WIN32
	pDstSeq = (UINT8 *) realloc(m_pSeqs2Assemb,memreq);
#else
	pDstSeq = (UINT8 *)mremap(m_pSeqs2Assemb,m_AllocMemSeqs2Assemb,memreq,MREMAP_MAYMOVE);
	if(pDstSeq == MAP_FAILED)
		pDstSeq = NULL;
#endif
	if(pDstSeq == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddSeq: Memory re-allocation to %lld bytes - %s",memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_AllocMemSeqs2Assemb = memreq;
	m_pSeqs2Assemb = pDstSeq;
	}

pDstSeq = &m_pSeqs2Assemb[m_Seqs2AssembLen];
*(UINT32 *)pDstSeq = SeqLen;
pDstSeq += sizeof(UINT32);
for(Idx = 0; Idx < SeqLen; Idx++)
	*pDstSeq++ = *pSeq++ & 0x07;
m_Seqs2AssembLen += SeqLen + sizeof(UINT32);

m_NumSeqs2Assemb += 1;
m_TotSeqs2AssembBases += SeqLen;
ReleaseSerialise();
return(eBSFSuccess);
}

// generate sfx over current concatenated assembled contigs
teBSFrsltCodes
CHomozyReduce::GenRdsSfx(void)
{
int ElSize;
UINT64 MaxSuffixEls;
UINT32 *pTmp;
MaxSuffixEls = m_RdsSfxHdr.ConcatSeqLen;
ElSize = sizeof(UINT32);
if(MaxSuffixEls > cMaxSfxBlkEls)
	ElSize += 1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Now generating suffix array with %lld index elements size %d bytes..",MaxSuffixEls,ElSize);

size_t ReqAllocMem;
ReqAllocMem = (size_t)((MaxSuffixEls + 10) * (INT64)ElSize);
if(m_pSuffixArray == NULL || m_AllocMemSfx == 0)
	{
	m_AllocMemSfx = ReqAllocMem;
#ifdef _WIN32
	m_pSuffixArray = (UINT32 *) malloc((size_t)m_AllocMemSfx);
	if(m_pSuffixArray == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Suffix array memory allocation of %lld bytes - %s",(INT64)m_AllocMemSfx,strerror(errno));
		m_AllocMemSfx = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pSuffixArray = (UINT32 *)mmap(NULL,m_AllocMemSfx, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pSuffixArray == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Suffix array memory allocation of %lld bytes through mmap()  failed - %s",(INT64)m_AllocMemSfx,strerror(errno));
		m_pSuffixArray = NULL;
		m_AllocMemSfx = 0;
		Reset(false);
		return(eBSFerrMem);
		}
#endif
	SetMaxMemWorkSetSize((size_t)(m_AllocMemConcat + m_AllocMemSeqs2Assemb + m_AllocMemSfx));
	}
else
	{

	// if worth the cost and suffix array can be reduced then do so - memory could be in short supply!
	if(ReqAllocMem >  (size_t)m_AllocMemSfx || ((ReqAllocMem * 10) <  ((size_t)m_AllocMemSfx * 11)))
		{
#ifdef _WIN32
		pTmp = (UINT32 *) realloc(m_pSuffixArray,ReqAllocMem);
#else
		pTmp = (UINT32 *)mremap(m_pSuffixArray,m_AllocMemSfx,ReqAllocMem,MREMAP_MAYMOVE);
		if(pTmp == MAP_FAILED)
			pTmp = NULL;
#endif
		if(pTmp == NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenRdsSfx: Memory re-allocation to %ld bytes - %s",ReqAllocMem,strerror(errno));
			return(eBSFerrMem);
			}
		m_AllocMemSfx = ReqAllocMem;
		m_pSuffixArray = pTmp;
		SetMaxMemWorkSetSize((size_t)(m_AllocMemConcat + m_AllocMemSeqs2Assemb + m_AllocMemSfx));
		}
	}


if(ElSize == 5)
	{
	UINT8 *pArr5 = (UINT8 *)m_pSuffixArray;
	for(UINT64 El = 0; El < MaxSuffixEls; El++)
		pArr5 = Pack5(El,pArr5);
	*pArr5 = cCSeqEOS;
	m_xpConcatSeqs = &m_pConcatSeqs[1];
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Now sorting...");
	m_MTqsort.qsort(m_pSuffixArray,(INT64)MaxSuffixEls,5,Sfx5SortFunc);
	}
else
	{
	for(UINT32 El = 0; El < MaxSuffixEls; El++)
		m_pSuffixArray[El] = El;
	m_pConcatSeqs[MaxSuffixEls] = cCSeqEOS;
	m_xpConcatSeqs = &m_pConcatSeqs[1];
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Now sorting...");
	m_MTqsort.qsort(m_pSuffixArray,(INT64)MaxSuffixEls,sizeof(UINT32),SfxSortFunc);
	}

m_NumSuffixEls = MaxSuffixEls;
m_SfxElSize = ElSize;


#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"GenRdsSfx: Suffix array generation completed");
return(eBSFSuccess);
}

// ChunkedWrite
// Seeks to specified 64bit file offset and writes to disk as chunks of no more than INT_MAX/16
teBSFrsltCodes
CHomozyReduce::ChunkedWrite(INT64 WrtOfs,UINT8 *pData,INT64 WrtLen)
{
return(ChunkedWrite(m_hOutFile,m_szOutFile,WrtOfs,pData,WrtLen));
}

teBSFrsltCodes
CHomozyReduce::ChunkedWrite(int hFile,char *pszFile,INT64 WrtOfs,UINT8 *pData,INT64 WrtLen)
{
int BlockLen;
if(_lseeki64(hFile,WrtOfs,SEEK_SET) != WrtOfs)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to seek to %ld on file %s - error %s",WrtOfs,pszFile,strerror(errno));
	return(eBSFerrFileAccess);
	}

while(WrtLen)
	{
	BlockLen = WrtLen > (INT64)(INT_MAX/16) ? (INT_MAX/16) : (int)WrtLen;
	WrtLen -= BlockLen;
	if(write(hFile,pData,BlockLen)!=BlockLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write to disk on file %s - error %s",pszFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
teBSFrsltCodes
CHomozyReduce::ChunkedRead(INT64 RdOfs,UINT8 *pData,INT64 RdLen)
{
return(ChunkedRead(m_hInFile,m_szInFile,RdOfs,pData,RdLen));
}

// ChunkedRead
// Seeks to specified 64bit file offset and reads from disk as chunks of no more than INT_MAX/32
teBSFrsltCodes
CHomozyReduce::ChunkedRead(int hFile,char *pszFile,INT64 RdOfs,UINT8 *pData,INT64 RdLen)
{
int BlockLen;
if(_lseeki64(hFile,RdOfs,SEEK_SET) != RdOfs)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to seek to %ld on file - error %s",RdOfs,pszFile,strerror(errno));
	return(eBSFerrFileAccess);
	}

while(RdLen)
	{
	BlockLen = RdLen > (INT64)(INT_MAX/32) ? (INT_MAX/32) : (int)RdLen;
	RdLen -= BlockLen;
	if(read(hFile,pData,BlockLen)!=BlockLen)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to read from disk on file %s - error %s",pszFile,strerror(errno));
		return(eBSFerrFileAccess);
		}
	pData += BlockLen;
	}
return(eBSFSuccess);
}


// SortContigsByLen
// Sort contigs by length
teBSFrsltCodes
CHomozyReduce::SortContigsByLen(void)
{
UINT32 ReadSeqID;
tsSfxdSeq *pRead;
UINT8 *pSeq;
UINT64 XFormID;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Sorting %d sequences ready for duplicate detection ...",m_NumSfxdSeqs);

	// first sort the reads
m_xpConcatSeqs = &m_pConcatSeqs[1];
m_MTqsort.qsort(m_pSfxdSeqs,(INT64)m_NumSfxdSeqs,sizeof(tsSfxdSeq),ContigsSortFunc);

	// next assign ReadSeqID's
UINT64 TotReadLens = 0;
pRead = m_pSfxdSeqs;
for(ReadSeqID = 1; ReadSeqID <= m_NumSfxdSeqs; ReadSeqID += 1, pRead += 1)
	{
	pRead->SeqID = ReadSeqID;
	TotReadLens += (UINT64)pRead->ReadLen;
	pSeq = &m_pConcatSeqs[pRead->ConcatSeqOfs];
	XFormID = IDtoXForm(pRead->SeqID);
	*(UINT32 *)pSeq = (UINT32)XFormID;
	pSeq += sizeof(UINT32);
	*pSeq = (UINT8)(XFormID >> 32);
	}

m_MeanReadLen = (int)(TotReadLens / m_NumSfxdSeqs);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed sorting %d sequences with mean length of %d",m_NumSfxdSeqs, m_MeanReadLen);
return(eBSFSuccess);
}



int  // SfxSortFunc for UINT32 suffix elements
CHomozyReduce::SfxSortFunc(const void *arg1, const void *arg2)
{
etSeqBase Base1;
etSeqBase Base2;
UINT8 Byte1;
UINT8 Byte2;
int MaxCmpLen = cAbsMaxCoreLen * 500;

UINT8 *pSeq1;
UINT8 *pSeq2;
UINT32 SeqIdx1 = *(UINT32 *)arg1;
UINT32 SeqIdx2 = *(UINT32 *)arg2;
pSeq1 = &m_xpConcatSeqs[SeqIdx1];
pSeq2 = &m_xpConcatSeqs[SeqIdx2];
do {
	Base1 = (Byte1 = *pSeq1++) & 0x07;
	Base2 = (Byte2 = *pSeq2++) & 0x07;
	if(Byte1 >= cCSeqSep)
		Base1 = Byte1;
	if(Byte2 >= cCSeqSep)
		Base2 = Byte2;
	if(Base1 < Base2)
		return(-1);
	if(Base1 > Base2)
		return(1);
	}
while(MaxCmpLen-- && (Base1 < 0x07 && Base2 < 0x07));
return(0);
}


int  // SfxSortFunc for 5byte suffix elements
CHomozyReduce::Sfx5SortFunc(const void *arg1, const void *arg2)
{
etSeqBase Base1;
etSeqBase Base2;
UINT8 Byte1;
UINT8 Byte2;
int MaxCmpLen = cAbsMaxCoreLen * 500;

UINT8 *pSeq1;
UINT8 *pSeq2;
UINT64 SeqIdx1 = Unpack5((UINT8 *)arg1);
UINT64 SeqIdx2 = Unpack5((UINT8 *)arg2);
pSeq1 = &m_xpConcatSeqs[SeqIdx1];
pSeq2 = &m_xpConcatSeqs[SeqIdx2];
do {
	Base1 = (Byte1 = *pSeq1++) & 0x07;
	Base2 = (Byte2 = *pSeq2++) & 0x07;
	if(Byte1 >= cCSeqSep)
		Base1 = Byte1;
	if(Byte2 >= cCSeqSep)
		Base2 = Byte2;
	if(Base1 < Base2)
		return(-1);
	if(Base1 > Base2)
		return(1);
	}
while(MaxCmpLen-- && (Base1 < 0x07 && Base2 < 0x07));
return(0);
}


// ContigsSortFunc
// Sort function for sorting reads by sequence
int
CHomozyReduce::ContigsSortFunc(const void *arg1, const void *arg2)
{
tsSfxdSeq *pR1 = (tsSfxdSeq *)arg1;
tsSfxdSeq *pR2 = (tsSfxdSeq *)arg2;
etSeqBase Base1;
etSeqBase Base2;
UINT8 Byte1;
UINT8 Byte2;

UINT8 *pSeq1;
UINT8 *pSeq2;
UINT64 SeqIdx1 = pR1->ConcatSeqOfs;
UINT64 SeqIdx2 = pR2->ConcatSeqOfs;
pSeq1 = &m_xpConcatSeqs[SeqIdx1 + 5];			// pts to 1st base of sequence
pSeq2 = &m_xpConcatSeqs[SeqIdx2 + 5];			// pts to 1st base of sequence
do {
	Base1 = (Byte1 = *pSeq1++) & 0x07;
	Base2 = (Byte2 = *pSeq2++) & 0x07;
	if(Byte1 >= cCSeqSep)
		Base1 = Byte1;
	if(Byte2 >= cCSeqSep)
		Base2 = Byte2;
	if(Base1 < Base2)
		return(-1);
	if(Base1 > Base2)
		return(1);
	}
while(Base1 < 0x07 && Base2 < 0x07);
return(0);
}