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


CProcRawReads::CProcRawReads(void)
{
m_hInFile = -1;
m_hOutFile = -1;
m_pWrtBuff = NULL;
m_pRdsBuff = NULL;
m_pDimerCnts = NULL;
m_pTrimerCnts = NULL;
m_pTetramerCnts = NULL;
m_pDistQualScores = NULL;
m_ppReadsIdx = NULL;
m_pDataBuff = NULL;
m_bActivity = false;
m_MTqsort.SetMaxThreads(cDfltSortThreads);
Reset();
}


CProcRawReads::~CProcRawReads(void)
{
Reset();
}

void
CProcRawReads::Reset(void)
{
if(m_hInFile != -1)
	{
	close(m_hInFile);
	m_hInFile = -1;
	}
if(m_hOutFile != -1)
	{
	close(m_hOutFile);
	m_hOutFile = -1;
	}

if(m_pWrtBuff != NULL)
	{
	delete m_pWrtBuff;
	m_pWrtBuff = NULL;
	}

if(m_pRdsBuff != NULL)
	{
	delete m_pRdsBuff;
	m_pRdsBuff = NULL;
	}
if(m_pDimerCnts != NULL)
	{
	delete m_pDimerCnts;
	m_pDimerCnts = NULL;
	}

if(m_pTrimerCnts != NULL)
	{
	delete m_pTrimerCnts;
	m_pTrimerCnts = NULL;
	}

if(m_pTetramerCnts != NULL)
	{
	delete m_pTetramerCnts;
	m_pTetramerCnts = NULL;
	}

if(m_pDistQualScores != NULL)
	{
	delete m_pDistQualScores;
	m_pDistQualScores = NULL;
	}

if(m_ppReadsIdx != NULL)
	{
#ifdef _WIN32
	free(m_ppReadsIdx);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_ppReadsIdx != MAP_FAILED)
		munmap(m_ppReadsIdx,m_ReadsIdxBuffAllocMem);
#endif
	m_ppReadsIdx = NULL;
	}


if(m_pDataBuff != NULL)
	{
#ifdef _WIN32
	free(m_pDataBuff);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pDataBuff != MAP_FAILED)
		munmap(m_pDataBuff,m_DataBuffAllocMem);
#endif
	m_pDataBuff = NULL;
	}
m_NumDescrReads = 0;
m_DataBuffAllocMem = 0;
m_DataBuffAllocMem = 0;
m_DataBuffOfs = 0;
m_bIsSOLiD = false;

memset(&m_FileHdr,0,sizeof(m_FileHdr));
}

void
CProcRawReads::ReportActivity(bool bActivity)
{
m_bActivity = bActivity;
}

// Hdr2Disk
// Writes file header out to disk
teBSFrsltCodes
CProcRawReads::Hdr2Disk(char *pszRdsFile)
{
int WrtLen;

m_FileHdr.Magic[0] = 'b'; m_FileHdr.Magic[1] = 'i'; m_FileHdr.Magic[2] = 'o'; m_FileHdr.Magic[3] = 'r';
m_FileHdr.Version = cBSFRRRdsVersion;
WrtLen = sizeof(tsBSFRdsHdr);
m_FileHdr.RdsOfs = WrtLen;

if(_lseeki64(m_hOutFile,0,SEEK_SET) ||
		!CUtility::SafeWrite(m_hOutFile,&m_FileHdr,WrtLen))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to write file header to disk file '%s'  - error %s",pszRdsFile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
return(eBSFSuccess);
}

teBSFrsltCodes
CProcRawReads::Disk2Hdr(char *pszRdsFile)
{
if(_lseeki64(m_hInFile,0,SEEK_SET)!=0)			// read in header..
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Seek failed to offset 0 - %s",pszRdsFile,strerror(errno));
	Reset();			// closes opened file..
	return(eBSFerrFileAccess);
	}

if(sizeof(tsBSFRdsHdr) != read(m_hInFile,&m_FileHdr,sizeof(tsBSFRdsHdr)))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Read of file header failed on %s - %s",pszRdsFile,strerror(errno));
	Reset();			// closes opened file..
	return(eBSFerrFileAccess);
	}

// header read, validate it as being a reads file header
if(tolower(m_FileHdr.Magic[0]) != 'b' ||
	tolower(m_FileHdr.Magic[1]) != 'i' ||
	tolower(m_FileHdr.Magic[2]) != 'o' ||
	tolower(m_FileHdr.Magic[3]) != 'r')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"%s opened but no magic signature - not a reads file",pszRdsFile);
	Reset();			// closes opened file..
	return(eBSFerrNotBioseq);
	}

	// can we handle this version?
if(m_FileHdr.Version < cBSFRRRdsVersionBack || m_FileHdr.Version > cBSFRRRdsVersion)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"%s opened as a pre-processed reads file - expected between version %d and %d, file version is %d",
					pszRdsFile,cBSFRRRdsVersionBack,cBSFRRRdsVersion,m_FileHdr.Version);
	Reset();			// closes opened file..
	return(eBSFerrFileVer);
	}

return(eBSFSuccess);
}

teBSFrsltCodes
CProcRawReads::GenDump(etPRRMode PMode,						// processing mode
		int NumReadsLimit,					// limit processing to this many reads
		char *pszInfile,					// input file
		char *pszOutfile)					// output into this file
{
int Rslt;
int BuffLen;
int BuffOfs;
int NumReadsProc;						// number of reads processed
int MaxLengthRead;						// max length read thus  far processed
int RdLen;								// last read() result when filling m_pRdsBuff[]
UINT8 *pSeqVal;
etSeqBase *pSeqFwd;
tsRawReadV5 *pReadV5;						// preprocessed reads could be either V5
tsRawReadV6 *pReadV6;						//  or V6

etSeqBase Sequence[cMaxFastQSeqLen];	// to hold sequence (sans quality scores) for current read
UINT8 Scores[cMaxFastQSeqLen];			// to hold quality scores for current read
UINT8 *pScore;
int SeqIdx;
unsigned long Elapsed;
unsigned long CurElapsed;

int SizeOfRawRead;
int CurReadLen;
int CurDescrLen;

CStopWatch CurTime;
CurTime.Start();
Elapsed = CurTime.ReadUSecs();

Reset();

if((m_pRdsBuff = new UINT8 [cRRRdsBuffAlloc])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes memory for reads buffering",cRRRdsBuffAlloc,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

m_pWrtBuff = new UINT8[cRRWrtBuffSize];
if(m_pWrtBuff == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Memory allocation of %d bytes for write buffer - %s",cRRWrtBuffSize,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

#ifdef _WIN32
m_hInFile = open(pszInfile, O_READSEQ ); // file access is normally sequential..
#else
m_hInFile = open64(pszInfile, O_READSEQ ); // file access is normally sequential..
#endif

if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszInfile,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

// expecting a preprocessed .rds file as input, header processing will confirm!
if((Rslt=Disk2Hdr(pszInfile))!=eBSFSuccess)
	{
	Reset();
	return((teBSFrsltCodes)Rslt);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reads file '%s' generator version: %d", pszInfile,m_FileHdr.Version);
if(m_FileHdr.Version < 4)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contains %d reads, PMode was %d, Quality was %d, %d bases 5' trimmed, %d bases 3' trimmed",
		m_FileHdr.NumRds,m_FileHdr.PMode,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);
	}
else	// version >= 4
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contains %d reads from %d source reads with duplicate sequences %s, PMode was %d, Quality was %d, %d bases 5' trimmed, %d bases 3' trimmed",
		m_FileHdr.NumRds,m_FileHdr.OrigNumReads,m_FileHdr.FlagsK ? "retained":"removed", m_FileHdr.PMode,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);

	}

if(m_FileHdr.FlagsCS == 1)
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reads were processed SOLiD colorspace files");
gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reads were processed from %d files",m_FileHdr.NumFiles);
char *pszSrcFile = (char *)m_FileHdr.FileNames;
for(BuffOfs=0;BuffOfs<m_FileHdr.NumFiles;BuffOfs++)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source file: '%s'",pszSrcFile);
	pszSrcFile += strlen(pszSrcFile) + 1;
	}


// ensure there are reads to generate Fastq or Fasta for
if(m_FileHdr.NumRds == 0 || m_FileHdr.TotReadsLen == 0)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Nothing to do, '%s' contains no reads",pszInfile);
	Reset();
	return(eBSFSuccess);
	}

// can now output iterated results
#ifdef _WIN32
m_hOutFile = open(pszOutfile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
if((m_hOutFile = open(pszOutfile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutfile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: unable to create/truncate output file '%s'",pszOutfile);
	Reset();
	return(eBSFerrCreateFile);
	}

int WrtOfs;
WrtOfs = 0;

// iterate each read sequence starting from the first
// tsRawRead struct size is version dependent
SizeOfRawRead = m_FileHdr.Version < 6 ? sizeof(tsRawReadV5) : sizeof(tsRawReadV6);

lseek(m_hInFile,(long)m_FileHdr.RdsOfs,SEEK_SET);
BuffLen = 0;
BuffOfs = 0;
NumReadsProc = 0;
MaxLengthRead = 0;
Elapsed = CurTime.ReadUSecs();
while((RdLen = read(m_hInFile,&m_pRdsBuff[BuffLen],cRRRdsBuffAlloc - BuffLen)) > 0)
	{
	BuffLen += RdLen;
	while((BuffLen - BuffOfs) >=  SizeOfRawRead)
		{
		if(m_FileHdr.Version < 6)
			{
			pReadV5 = (tsRawReadV5 *)&m_pRdsBuff[BuffOfs];
			pSeqVal = &pReadV5->Read[pReadV5->DescrLen+1];
			CurReadLen = pReadV5->ReadLen;
			CurDescrLen = pReadV5->DescrLen;
			}
		else
			{
			pReadV6 = (tsRawReadV6 *)&m_pRdsBuff[BuffOfs];
			pSeqVal = &pReadV6->Read[pReadV6->DescrLen+1];
			CurReadLen = pReadV6->ReadLen;
			CurDescrLen = pReadV6->DescrLen;
			}

		if((int)(CurReadLen + CurDescrLen + SizeOfRawRead) > (BuffLen - BuffOfs))
				break;
		BuffOfs += CurReadLen + CurDescrLen + SizeOfRawRead;
		if(m_bActivity && !(NumReadsProc % 1000))		// reduces the overhead of checking on elapsed time if it is done every 1000 reads...
			{
			CurElapsed = CurTime.ReadUSecs();
			if((CurElapsed - Elapsed) > 60)
				{
				Elapsed = CurElapsed;
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing read %d",NumReadsProc);
				}
			}


		if(MaxLengthRead < (int)CurReadLen)
			MaxLengthRead = CurReadLen;

		// get descriptor, sequence and quality values
		pSeqFwd = Sequence;
		pScore = Scores;
		for(SeqIdx = 0; SeqIdx < CurReadLen; SeqIdx++)
			{
			*pSeqFwd++ = *pSeqVal & 0x07;
			*pScore++ = (*pSeqVal++ >> 3) & 0x01f;
			}

		if((WrtOfs + (cRRWrtBuffSize/8)) > cRRWrtBuffSize)
			{
			if(!CUtility::SafeWrite(m_hOutFile,m_pWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenDump: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}

		if(m_FileHdr.Version < 6)
			WrtOfs += sprintf((char *)&m_pWrtBuff[WrtOfs],">lcl|%s %d|%s|%d|%d|%d|%d|%d|%d\n",
						   (char *)pReadV5->Read,pReadV5->ReadID, m_FileHdr.PMode == ePMRRNewPaired ? "p" : "s",
							pReadV5->FileID,pReadV5->ReadLen,pReadV5->NumReads,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);
		else
			WrtOfs += sprintf((char *)&m_pWrtBuff[WrtOfs],">lcl|%s %d|%s|%d|%d|%d|%d|%d|%d\n",
						   (char *)pReadV6->Read,pReadV6->ReadID, m_FileHdr.PMode == ePMRRNewPaired ? "p" : "s",
							pReadV6->FileID,pReadV6->ReadLen,pReadV6->NumReads,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);


		int NumCols;
		int SeqOfs = 0;
		int Len = CurReadLen;
		while(Len)
			{
			NumCols = Len > 70 ? 70 : Len;
			CSeqTrans::MapSeq2Ascii(&Sequence[SeqOfs],NumCols,(char *)&m_pWrtBuff[WrtOfs]);
			WrtOfs += NumCols;
			WrtOfs += sprintf((char *)&m_pWrtBuff[WrtOfs],"\n");
			Len -= NumCols;
			SeqOfs += NumCols;
			}
		if((WrtOfs + cRRWrtBuffSize/8) > cRRWrtBuffSize)
			{
			if(!CUtility::SafeWrite(m_hOutFile,m_pWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenDump: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}
		NumReadsProc+=1;
		if(NumReadsLimit && NumReadsProc >= NumReadsLimit)
			break;
		}
	if(NumReadsLimit && NumReadsProc >= NumReadsLimit)
		break;
	memmove(m_pRdsBuff,&m_pRdsBuff[BuffOfs],BuffLen - BuffOfs);
	BuffLen -= BuffOfs;
	BuffOfs = 0;
	}
if(WrtOfs)
	{
	if(!CUtility::SafeWrite(m_hOutFile,m_pWrtBuff,WrtOfs))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenDump: Write to %s - %s",pszOutfile,strerror(errno));
		Reset();
		return(eBSFerrFileAccess);
		}
	}
Reset();
return(eBSFSuccess);
}

teBSFrsltCodes
CProcRawReads::GenStats(etPRRMode PMode,	// processing mode
		char *pszInfile,					// input file
		char *pszOutfile)					// output into this file
{
int Rslt;
int BuffLen;
int BuffOfs;
int NumReadsProc;						// number of reads processed
int MaxLengthRead;						// max length read thus  far processed
int RdLen;								// last read() result when filling m_pRdsBuff[]
UINT8 *pSeqVal;
etSeqBase *pSeqFwd;
tsRawReadV5 *pReadV5;						// current preprocessed read being processed if V5
tsRawReadV6 *pReadV6;						// current preprocessed read being processed if V6
etSeqBase Sequence[cMaxFastQSeqLen];	// to hold sequence (sans quality scores) for current read
UINT8 Scores[cMaxFastQSeqLen];			// to hold quality scores for current read
UINT8 *pScore;
int SeqIdx;

UINT32 BaseCnts[6];						// total counts for each base including 'N's
UINT32 DistBaseCnts[6][cMaxFastQSeqLen]; // distribution of base counts along reads
UINT32 SeqNsCnts[cMaxFastQSeqLen];		// distribution of number of 'N's over all reads
int CntIdx;
int SeqNsIdx;

int SizeOfRawRead;
int CurReadLen;
int CurDescrLen;

unsigned long Elapsed;
unsigned long CurElapsed;
CStopWatch CurTime;
CurTime.Start();
Elapsed = CurTime.ReadUSecs();

Reset();

if((m_pDistQualScores = new UINT32 [cMaxFastQSeqLen * 2 * 0x01f])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for quality counts",strerror(errno));
	return(eBSFerrMem);
	}

if((m_pDimerCnts = new UINT32 [5*5*cMaxFastQSeqLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for dimer counts",strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

if((m_pTrimerCnts = new UINT32 [5*5*5*cMaxFastQSeqLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for trimer counts",strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

if((m_pTetramerCnts = new UINT32 [5*5*5*5*cMaxFastQSeqLen])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for tetramer counts",strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

memset(BaseCnts,0,sizeof(BaseCnts));
memset(m_pDimerCnts,0,sizeof(UINT32) * 5*5*cMaxFastQSeqLen);
memset(m_pTrimerCnts,0,sizeof(UINT32) * 5*5*5*cMaxFastQSeqLen);
memset(m_pTetramerCnts,0,sizeof(UINT32) * 5*5*5*5*cMaxFastQSeqLen);


memset(SeqNsCnts,0,sizeof(SeqNsCnts));
memset(DistBaseCnts,0,sizeof(DistBaseCnts));
memset(m_pDistQualScores,0,sizeof(UINT32) * cMaxFastQSeqLen * 2 * 0x01f);

if((m_pRdsBuff = new UINT8 [cRRRdsBuffAlloc])==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate %d bytes memory for reads buffering",cRRRdsBuffAlloc,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

#ifdef _WIN32
m_hInFile = open(pszInfile, O_READSEQ ); // file access is normally sequential..
#else
m_hInFile = open64(pszInfile, O_READSEQ ); // file access is normally sequential..
#endif

if(m_hInFile == -1)							// check if file open succeeded
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to open %s - %s",pszInfile,strerror(errno));
	Reset();
	return(eBSFerrOpnFile);
	}

// expecting a preprocessed .rds file as input, header processing will confirm!
if((Rslt=Disk2Hdr(pszInfile))!=eBSFSuccess)
	{
	Reset();
	return((teBSFrsltCodes)Rslt);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reads file '%s' generator version: %d", pszInfile,m_FileHdr.Version);
if(m_FileHdr.FlagsCS == 1)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Reads were processed from SOLiD colorspace");

if(m_FileHdr.Version < 4)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contains %d reads, PMode was %d, Quality was %d, %d bases 5' trimmed, %d bases 3' trimmed",
		m_FileHdr.NumRds,m_FileHdr.PMode,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);
	}
else
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Contains %d reads from %d source reads with duplicate sequences %s, PMode was %d, Quality was %d, %d bases 5' trimmed, %d bases 3' trimmed",
		m_FileHdr.NumRds,m_FileHdr.OrigNumReads,m_FileHdr.FlagsK ? "retained":"removed", m_FileHdr.PMode,m_FileHdr.QMode,m_FileHdr.Trim5,m_FileHdr.Trim3);

	}

gDiagnostics.DiagOutMsgOnly(eDLInfo,"Reads were processed from %d files",m_FileHdr.NumFiles);
char *pszSrcFile = (char *)m_FileHdr.FileNames;
for(BuffOfs=0;BuffOfs<m_FileHdr.NumFiles;BuffOfs++)
	{
	gDiagnostics.DiagOutMsgOnly(eDLInfo,"Source file: '%s'",pszSrcFile);
	pszSrcFile += strlen(pszSrcFile) + 1;
	}


// ensure there are reads to generate stats for
if(m_FileHdr.NumRds == 0 || m_FileHdr.TotReadsLen == 0)
	{
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Nothing to do, '%s' contains no reads",pszInfile);
	Reset();
	return(eBSFSuccess);
	}

	// iterate each read sequence starting from the first
// tsRawRead struct size is version dependent
SizeOfRawRead = m_FileHdr.Version < 6 ? sizeof(tsRawReadV5) : sizeof(tsRawReadV6);

lseek(m_hInFile,(long)m_FileHdr.RdsOfs,SEEK_SET);
BuffLen = 0;
BuffOfs = 0;
NumReadsProc = 0;
MaxLengthRead = 0;
Elapsed = CurTime.ReadUSecs();
while((RdLen = read(m_hInFile,&m_pRdsBuff[BuffLen],cRRRdsBuffAlloc - BuffLen)) > 0)
	{
	BuffLen += RdLen;
	while((BuffLen - BuffOfs) >=  sizeof(tsRawReadV5))
		{
		if(m_FileHdr.Version < 6)
			{
			pReadV5 = (tsRawReadV5 *)&m_pRdsBuff[BuffOfs];
			pSeqVal = &pReadV5->Read[pReadV5->DescrLen+1];
			CurReadLen = pReadV5->ReadLen;
			CurDescrLen = pReadV5->DescrLen;
			}
		else
			{
			pReadV6 = (tsRawReadV6 *)&m_pRdsBuff[BuffOfs];
			pSeqVal = &pReadV6->Read[pReadV6->DescrLen+1];
			CurReadLen = pReadV6->ReadLen;
			CurDescrLen = pReadV6->DescrLen;
			}

		if((int)(CurReadLen + CurDescrLen + SizeOfRawRead) > (BuffLen - BuffOfs))
				break;
		BuffOfs += CurReadLen + CurDescrLen + SizeOfRawRead;

		if(m_bActivity && !(NumReadsProc % 1000))		// reduces the overhead of checking on elapsed time if it is done every 1000 reads...
			{
			CurElapsed = CurTime.ReadUSecs();
			if((CurElapsed - Elapsed) > 60)
				{
				Elapsed = CurElapsed;
				gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing read %d",NumReadsProc);
				}
			}

		if(MaxLengthRead < (int)CurReadLen)
			MaxLengthRead = CurReadLen;

		// get separate sequence and quality values
		pSeqFwd = Sequence;
		pScore = Scores;
		for(SeqIdx = 0; SeqIdx < CurReadLen; SeqIdx++)
			{
			*pSeqFwd++ = *pSeqVal & 0x07;
			*pScore++ = (*pSeqVal++ >> 3) & 0x01f;
			}

		SeqNsIdx = 0;
		CntIdx  = 0;
		pSeqFwd = Sequence;
		pScore = Scores;
		for(SeqIdx = 0; SeqIdx < CurReadLen; SeqIdx++,pSeqFwd++,pScore++)
			{
			if(*pSeqFwd == eBaseN)
				SeqNsIdx += 1;

			int Offs = (*pScore * cMaxFastQSeqLen) + SeqIdx;
			m_pDistQualScores[Offs] += 1;

			BaseCnts[*pSeqFwd] += 1;
			DistBaseCnts[*pSeqFwd][SeqIdx] += 1;

			if(SeqIdx < CurReadLen-1)
				{
				CntIdx = (((*pSeqFwd * 5) + pSeqFwd[1]) * cMaxFastQSeqLen)+SeqIdx;
				m_pDimerCnts[CntIdx] += 1;
				}

			if(SeqIdx < CurReadLen-2)
				{
				CntIdx = (((((*pSeqFwd * 5) + pSeqFwd[1])*5) + pSeqFwd[2]) * cMaxFastQSeqLen)+SeqIdx;
				m_pTrimerCnts[CntIdx] += 1;
				}

			if(SeqIdx < CurReadLen-3)
				{
				CntIdx = (((((((*pSeqFwd * 5) + pSeqFwd[1])*5) + pSeqFwd[2]) * 5) +  pSeqFwd[3]) * cMaxFastQSeqLen)+SeqIdx;
				m_pTetramerCnts[CntIdx] += 1;
				}
			}
		SeqNsCnts[SeqNsIdx] += 1;

#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif

		NumReadsProc+=1;
		}
	memmove(m_pRdsBuff,&m_pRdsBuff[BuffOfs],BuffLen - BuffOfs);
	BuffLen -= BuffOfs;
	BuffOfs = 0;
	}

close(m_hInFile);
m_hInFile = -1;

// confirm expected number of reads were processed
if(m_FileHdr.NumRds != NumReadsProc)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Expected to process %d reads, but processed %d from file '%s'",
								m_FileHdr.NumRds,NumReadsProc,pszOutfile);
	Reset();
	return(eBSFerrCreateFile);
	}


// can now output results
#ifdef _WIN32
m_hOutFile = open(pszOutfile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
if((m_hOutFile = open(pszOutfile,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutfile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: unable to create/truncate output file '%s'",pszOutfile);
	Reset();
	return(eBSFerrCreateFile);
	}


int ReadIdx;
int WrtOfs;
char szWrtBuff[cMaxReadLen+1];
WrtOfs = 0;

// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= MaxLengthRead; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	}
for(ReadIdx = 0; ReadIdx < 5; ReadIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(ReadIdx));
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < MaxLengthRead; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",DistBaseCnts[ReadIdx][SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
WrtOfs = 0;

// repeat monomer distribution but this time as proportions
for(SeqIdx = 1; SeqIdx <= MaxLengthRead; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	}
for(ReadIdx = 0; ReadIdx < 4; ReadIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(ReadIdx));
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < MaxLengthRead; SeqIdx++)
		{
		int ColTot = 0;
		for(int RowIdx = 0; RowIdx < 4; RowIdx++)
			ColTot += DistBaseCnts[RowIdx][SeqIdx];
		if(ColTot > 0)
			WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%1.4f",(double)DistBaseCnts[ReadIdx][SeqIdx]/(double)ColTot);
		else
			WrtOfs += sprintf(&szWrtBuff[WrtOfs],"0.0000");
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
WrtOfs = 0;

// dimer distributions
// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= MaxLengthRead-1; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	}

int Base1Idx;
for(Base1Idx = 0; Base1Idx < 25; Base1Idx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(Base1Idx / 5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(Base1Idx % 5));
	CntIdx = Base1Idx * cMaxFastQSeqLen;
	for(SeqIdx = 0; SeqIdx < MaxLengthRead-1; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_pDimerCnts[CntIdx + SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
WrtOfs = 0;

// trimer distributions
// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= MaxLengthRead-2; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	}


for(Base1Idx = 0; Base1Idx < (5*5*5); Base1Idx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii((Base1Idx/5)/5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii((Base1Idx/5)%5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(Base1Idx%5));
	CntIdx = Base1Idx * cMaxFastQSeqLen;
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < MaxLengthRead-2; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_pTrimerCnts[CntIdx+SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
WrtOfs = 0;

// tetramer distributions
// row 1 is the base position (1..MaxLengthRead) heading
for(SeqIdx = 1; SeqIdx <= MaxLengthRead-3; SeqIdx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",SeqIdx);
	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	}


for(Base1Idx = 0; Base1Idx < (5*5*5*5); Base1Idx++)
	{
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%c",CSeqTrans::MapBase2Ascii(((Base1Idx/5)/5)/5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(((Base1Idx/5)/5)%5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii((Base1Idx/5)%5));
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"%c",CSeqTrans::MapBase2Ascii(Base1Idx%5));
	CntIdx = Base1Idx * cMaxFastQSeqLen;

	if(WrtOfs + 100 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < MaxLengthRead-3; SeqIdx++)
		{
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_pTetramerCnts[CntIdx+SeqIdx]);
		if(WrtOfs + 100 >= sizeof(szWrtBuff))
			{
			if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}
		}
	}
WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n\n");
if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
WrtOfs = 0;

// output the quality scores distribution over the read length
// subsequent rows
// scale phred back to 0..40
int Phred;
// if wanting to convert to the original probability then use
//Prob = 1.0/(pow(10.0,(ReadIdx * 40.0/31.0)/10.0))
for(ReadIdx = 0; ReadIdx <= 0x01f; ReadIdx++)
	{
	//Prob = 1.0/(pow(10.0,(ReadIdx * 40.0/31.0)/10.0));
	Phred = (ReadIdx * 40)/31;
	WrtOfs += sprintf(&szWrtBuff[WrtOfs],"\n%2.2d",Phred);
	if(WrtOfs + 16 >= sizeof(szWrtBuff))
		{
		if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	for(SeqIdx = 0; SeqIdx < MaxLengthRead; SeqIdx++)
		{
		int Offs = (ReadIdx * cMaxFastQSeqLen) + SeqIdx;
		WrtOfs += sprintf(&szWrtBuff[WrtOfs],",%d",m_pDistQualScores[Offs]);
		if(WrtOfs + 16 >= sizeof(szWrtBuff))
			{
			if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
				Reset();
				return(eBSFerrFileAccess);
				}
			WrtOfs = 0;
			}
		}
	}


if(WrtOfs)
	{
	if(!CUtility::SafeWrite(m_hOutFile,szWrtBuff,WrtOfs))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"GenStats: Write to %s - %s",pszOutfile,strerror(errno));
		Reset();
		return(eBSFerrFileAccess);
		}
	}

Reset(); // will close file handles..

return(eBSFSuccess);
}

teBSFrsltCodes
CProcRawReads::WriteToFile(etPRRMode PMode,
			etFQMethod Quality,					// fastq quality value method
			bool bKeepDups,
			int Trim5,
			int Trim3,
			char *pszOutFile)
{
int Rslt;
int Idx;
int WrtOfs;
UINT64 TotReadsLen;
tsRawReadV6 *pReadV6;

unsigned long Elapsed;
unsigned long CurElapsed;
CStopWatch CurTime;
CurTime.Start();
Elapsed = CurTime.ReadUSecs();

if(pszOutFile == NULL || pszOutFile[0] == '\0')
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"WriteToFile: no output file specified");
	return(eBSFerrCreateFile);
	}



#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
if((m_hOutFile = open(pszOutFile,O_WRONLY | O_CREAT, S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Starting to write to file '%s'",pszOutFile);

m_FileHdr.NumRds = 0;
m_FileHdr.TotReadsLen = 0;

// write out a small header
if((Rslt=Hdr2Disk(pszOutFile))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to write output file header '%s'",pszOutFile);
	Reset();
	return(eBSFerrCreateFile);
	}

m_pWrtBuff = new UINT8[cRRWrtBuffSize];
if(m_pWrtBuff == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Memory allocation of %d bytes for write buffer - %s",cRRWrtBuffSize,strerror(errno));
	Reset();
	return(eBSFerrMem);
	}

if(_lseeki64(m_hOutFile,m_FileHdr.RdsOfs,SEEK_SET)!=m_FileHdr.RdsOfs)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to seek to start of reads on disk file '%s'  - error %s",pszOutFile,strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}
WrtOfs = 0;
TotReadsLen = 0;
Elapsed = CurTime.ReadUSecs();
for(Idx = 0; Idx < (int)m_NumDescrReads; Idx++)
	{
	pReadV6 = m_ppReadsIdx[Idx];
	if(pReadV6->ReadID == 0)
		continue;
	if((WrtOfs + 2 * (cMaxFastQSeqLen + cMaxDescrIDLen + sizeof(tsRawReadV6))) >= cRRWrtBuffSize)		// allow a small safety margin
		{
		if(!CUtility::SafeWrite(m_hOutFile,m_pWrtBuff,WrtOfs))
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",WrtOfs, pszOutFile, strerror(errno));
			Reset();
			return(eBSFerrFileAccess);
			}
		WrtOfs = 0;
		}
	m_FileHdr.NumRds += 1;
	pReadV6->ReadID =m_FileHdr.NumRds;
	memmove(&m_pWrtBuff[WrtOfs],pReadV6,sizeof(tsRawReadV6) + pReadV6->DescrLen + pReadV6->ReadLen);
	WrtOfs += sizeof(tsRawReadV6) + pReadV6->DescrLen + pReadV6->ReadLen;
	TotReadsLen += pReadV6->DescrLen + pReadV6->ReadLen + 1;

	if(!(m_FileHdr.NumRds % 1000))		// reduces the overhead of checking on elapsed time if it is done every 1000 reads...
		{
		CurElapsed = CurTime.ReadUSecs();
		if((CurElapsed - Elapsed) > 60)
			{
			Elapsed = CurElapsed;
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Written %d reads to file '%s'",m_FileHdr.NumRds,pszOutFile);
			}
		}
	}

if(WrtOfs && !CUtility::SafeWrite(m_hOutFile,m_pWrtBuff,WrtOfs))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors whilst writing %d bytes - '%s' - %s",WrtOfs, pszOutFile, strerror(errno));
	Reset();
	return(eBSFerrFileAccess);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Written %d reads to file '%s'",m_FileHdr.NumRds,pszOutFile);

// can now update header with actual numbers
m_FileHdr.OrigNumReads = m_NumDescrReads;
m_FileHdr.FlagsK = bKeepDups;
m_FileHdr.FlagsCS = m_bIsSOLiD;
m_FileHdr.FlagsPR = PMode == ePMRRNewPaired ? 1 : 0;
m_FileHdr.TotReadsLen = TotReadsLen;
m_FileHdr.PMode = (UINT8)PMode;
m_FileHdr.QMode = Quality;
m_FileHdr.Trim5 = Trim5;
m_FileHdr.Trim3 = Trim3;
if((Rslt=Hdr2Disk(pszOutFile))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to write output file header '%s'",pszOutFile);
	Reset();
	return(eBSFerrCreateFile);
	}

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Write to file '%s' completed",pszOutFile);

Reset();	// will close output file

return(eBSFSuccess);
}


teBSFrsltCodes
CProcRawReads::LoadAndProcessReads(etPRRMode PMode,		// processing mode
		UINT32 NumReadsLimit,				// limit processing to this many reads
		bool bKeepDups,						// true if duplicate reads not to be filtered out
		etFQMethod Quality,					// fastq quality value method
		int Trim5,							// trim this many bases off leading sequence 5' end
		int Trim3,							// trim this many bases off trailing sequence 3' end
		int NumInputFileSpecs,				// number of input file specs
		char *pszInfileSpecs[],				// names of inputs file (wildcards allowed unless in dump mode) containing raw reads
		char *pszInPairFile,				// if paired reads processing then file containing paired ends
		char *pszOutFile)					// output into this file
{
int Rslt;
char *pszInfile;
int Idx;			// general processing iteration index

Reset();

int NumInputFilesProcessed;

NumInputFilesProcessed = 0;
if(PMode != ePMRRNewPaired)
	{
	CSimpleGlob glob(SG_GLOB_FULLSORT);
	for(Idx = 0; Idx < NumInputFileSpecs; Idx++)
		{
		glob.Init();
		if(glob.Add(pszInfileSpecs[Idx]) < SG_SUCCESS)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInfileSpecs[Idx]);
			return(eBSFerrOpnFile);	// treat as though unable to open file
			}

		if(glob.FileCount() <= 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source raw reads file matching '%s",pszInfileSpecs[Idx]);
			continue;
			}

		Rslt = eBSFSuccess;
		for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
			{
			pszInfile = glob.File(FileID);
			NumInputFilesProcessed += 1;
			if(NumInputFilesProcessed > cRRMaxInFiles)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many input files (max allowed is %d)",cRRMaxInFiles);
				Reset();
				return(eBSFerrNumSrcFiles);
				}
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing reads from raw sequence file '%s'\n",pszInfile);
			Rslt = LoadReads(false,0,NumReadsLimit,Quality,Trim5,Trim3,NumInputFilesProcessed,pszInfile);
			if(Rslt != eBSFSuccess)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence file '%s'\n",pszInfile);
				Reset();
				return((teBSFrsltCodes)Rslt);
				}
			}
		}
	}
else
	{
	int ExpNumPairedReads;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading reads from paired end raw sequence file '%s'",pszInfileSpecs[0]);
	Rslt = LoadReads(false,1,NumReadsLimit,Quality,Trim5,Trim3,++NumInputFilesProcessed,pszInfileSpecs[0]);
	if(Rslt != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence paired reads file '%s'\n",pszInfileSpecs[0]);
		Reset();
		return((teBSFrsltCodes)Rslt);
		}
	ExpNumPairedReads = m_NumDescrReads;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Loading reads from partner paired end raw sequence file '%s'",pszInPairFile);
	Rslt = LoadReads(true,1,NumReadsLimit,Quality,Trim5,Trim3,++NumInputFilesProcessed,pszInPairFile);
	if(Rslt != eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence paired reads file '%s'\n",pszInPairFile);
		Reset();
		return((teBSFrsltCodes)Rslt);
		}
	// Ensure that all reads were paired...
	if(m_NumDescrReads != ExpNumPairedReads * 2)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: File '%s' contained %d reads, file '%s' contained %d reads - expected all reads to be paired",pszInfileSpecs[0],ExpNumPairedReads,pszInPairFile,m_NumDescrReads-ExpNumPairedReads);
		Reset();
		return((teBSFrsltCodes)Rslt);
		}

	}


#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d raw sequence files were accepted for filtering", NumInputFilesProcessed);
if(NumInputFilesProcessed == 0)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, no raw sequence files to be filtered");
	return(eBSFSuccess);
	}


// construct an index of ptrs to each and initialise
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Constructing index for %d reads\n",m_NumDescrReads);

if(m_ppReadsIdx == NULL)
	{
	size_t memreq = sizeof(tsRawReadV6 *) * m_NumDescrReads;
#ifdef _WIN32
	m_ppReadsIdx = (tsRawReadV6 **) malloc((size_t)memreq);
	if(m_ppReadsIdx == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Reads Idx: Memory allocation of %lld bytes failed",(INT64)memreq);
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_ppReadsIdx = (tsRawReadV6 **)mmap(NULL,(size_t)memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_ppReadsIdx == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %lld bytes through mmap()  failed",(INT64)memreq,strerror(errno));
		m_ppReadsIdx = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_ReadsIdxBuffAllocMem = memreq;
	}

tsRawReadV6 *pRead = (tsRawReadV6 *)m_pDataBuff;
tsRawReadV6 *pRead1;
UINT8 *pSeq;
size_t BuffOfs = 0;
for(Idx = 0; Idx < (int)m_NumDescrReads; Idx++)
	{
	m_ppReadsIdx[Idx] = pRead;
	BuffOfs += sizeof(tsRawReadV6) + pRead->ReadLen + pRead->DescrLen;
	pRead = (tsRawReadV6 *)&m_pDataBuff[BuffOfs];
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: sorting %d reads",m_NumDescrReads);

m_MTqsort.qsort(m_ppReadsIdx,(INT64)m_NumDescrReads,sizeof(tsRawReadV6 *),SortReads);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: sort of %d reads completed",m_NumDescrReads);

if(!bKeepDups)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: starting to identify those reads with unique sequences");

	int NumUniqueReadSeqs = 0;
	int Ofs;
	long SumQuals[cMaxFastQSeqLen];		// to hold sums of quality scores used to determine average if multiple reads have identical sequences
	pRead = m_ppReadsIdx[0];
	for(Idx = 1; Idx < (int)m_NumDescrReads; Idx++)
		{
		pRead1 = m_ppReadsIdx[Idx];
		Rslt = CompareRead(pRead,pRead1);
		if(Rslt > 0)							// problem if pRead ever compares as lexicographically higher than pRead1
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Unexpected read lexicographic order");
		if(Rslt != 0)
			{
			if(!bKeepDups && pRead->NumReads > 1)		// if there were multiple reads with identical sequences then use the averaged quality scores
				{
				pSeq = &pRead->Read[pRead->DescrLen+1];
				for(Ofs = 0; Ofs < pRead->ReadLen; Ofs++, pSeq++)
					{
					SumQuals[Ofs] += (*pSeq >> 3) & 0x01f;
					*pSeq = *pSeq & 0x07 | (UINT8)(((SumQuals[Ofs]+1) / pRead->NumReads) << 3);
					}
				}
			NumUniqueReadSeqs += 1;
			pRead = pRead1;
			continue;
			}

			// reads are identical, sum the quality scores so can later determine the average
		pSeq = &pRead1->Read[pRead->DescrLen+1];
		for(Ofs = 0; Ofs < pRead->ReadLen; Ofs++, pSeq++)
			if(pRead->NumReads == 1)
				SumQuals[Ofs] = (*pSeq >> 3) & 0x01f;
			else
				SumQuals[Ofs] += (*pSeq >> 3) & 0x01f;
		pRead1->ReadID = 0;		// this read no longer of relevance
		pRead->NumReads += 1;
		}

	if(pRead->NumReads > 1)		// if there were multiple reads with identical sequences then use averaged quality scores
		{
		pSeq = &pRead->Read[pRead->DescrLen+1];
		for(Ofs = 0; Ofs < pRead->ReadLen; Ofs++, pSeq++)
			{
			SumQuals[Ofs] += (*pSeq >> 3) & 0x01f;
			*pSeq = *pSeq & 0x07 | (UINT8)(((SumQuals[Ofs]+1) / pRead->NumReads) << 3);
			}
		pRead1->ReadID = 0;		// this read no longer of relevance
		}
	NumUniqueReadSeqs += 1;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: There are %d unique read sequences",NumUniqueReadSeqs);
	}

if(pszOutFile != NULL && pszOutFile[0] != '\0')
	{
	Rslt = WriteToFile(PMode,Quality,bKeepDups,Trim5,Trim3,pszOutFile);
	}
return((teBSFrsltCodes)Rslt);
}


// LocateFileName
// Locates and returns ptr to filename corresponding to FileID
char *
CProcRawReads::LocateFileName(int FileID)
{
int FileNameOfs;
char *pszFile;
if(m_FileHdr.NumFiles < FileID || FileID > m_FileHdr.NumFiles)
	return(NULL);

FileNameOfs = 0;
pszFile = (char *)&m_FileHdr.FileNames[0];
while(--FileID > 0)
	{
	FileNameOfs = (int)strlen(pszFile) + 1;
	pszFile += FileNameOfs;
	}
return(pszFile);
}



teBSFrsltCodes
CProcRawReads::LoadAndProcessReadsDE(etPRRMode PMode,						// processing mode
		UINT32 NumReadsLimit,				// limit processing to this many reads
		int Trim5,							// trim this many bases off leading sequence 5' end
		int Trim3,							// trim this many bases off trailing sequence 3' end
		int	MinSampleCnts,					// minimum sample counts
		int	MinTotalCnts,					// minimum total samples counts
		int NumInputFileSpecs,				// number of input file specs
		char *pszInfileSpecs[],				// names of inputs file (wildcards allowed unless in dump mode) containing raw reads
		char *pszOutFile)					// output into this file
{
int Rslt;
char *pszInfile;
int Idx;			// general processing iteration index
tsRawReadV6 *pRead;
tsRawReadV6 *pRead1;
size_t BuffOfs;
UINT32 SampleReadCnts[cRRMaxInFiles];		// read counts for each sample
char szBuff[(cMaxReadLen + 1) *4];
int NumUniqueReadSeqs;
int SeqIdx;
UINT8 *pSrc;
etSeqBase Sequence[cMaxReadLen+1];
etSeqBase *pDst = Sequence;
int SmplIdx;
bool bMinAnySampleCnts;
int NumSamplesTotalCnts;

Reset();

int NumInputFilesProcessed;

NumInputFilesProcessed = 0;

#ifdef _WIN32
m_hOutFile = open(pszOutFile,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE) );
#else
if((m_hOutFile = open(pszOutFile,O_WRONLY | O_CREAT, S_IREAD | S_IWRITE))!=-1)
    if(ftruncate(m_hOutFile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",pszOutFile,strerror(errno));
			return(eBSFerrCreateFile);
			}
#endif

if(m_hOutFile < 0)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output file '%s'",pszOutFile);
	Reset();
	return(eBSFerrCreateFile);
	}

CSimpleGlob glob(SG_GLOB_FULLSORT);
for(Idx = 0; Idx < NumInputFileSpecs; Idx++)
	{
	glob.Init();
	if(glob.Add(pszInfileSpecs[Idx]) < SG_SUCCESS)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to glob '%s",pszInfileSpecs[Idx]);
		return(eBSFerrOpnFile);	// treat as though unable to open file
		}

	if(glob.FileCount() <= 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to locate any source raw reads file matching '%s",pszInfileSpecs[Idx]);
		continue;
		}

	Rslt = eBSFSuccess;
	for (int FileID = 0; Rslt >= eBSFSuccess &&  FileID < glob.FileCount(); ++FileID)
		{
		pszInfile = glob.File(FileID);
		NumInputFilesProcessed += 1;
		if(NumInputFilesProcessed > cRRMaxInFiles)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Too many input files (max allowed is %d)",cRRMaxInFiles);
			Reset();
			return(eBSFerrNumSrcFiles);
			}

		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing reads from raw sequence file '%s'\n",pszInfile);
		Rslt = LoadReads(false,0,NumReadsLimit,eFQIgnore,Trim5,Trim3,NumInputFilesProcessed,pszInfile);
		if(Rslt != eBSFSuccess)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: Load failed for input raw sequence file '%s'\n",pszInfile);
			Reset();
			return((teBSFrsltCodes)Rslt);
			}
		}
	}



#ifdef _DEBUG
#ifdef _WIN32
_ASSERTE( _CrtCheckMemory());
#endif
#endif
gDiagnostics.DiagOut(eDLInfo,gszProcName,"%d raw sequence files were accepted for filtering", NumInputFilesProcessed);
if(NumInputFilesProcessed < 1)
	{
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Nothing to do, require at least one raw sequence files to be processed");
	return(eBSFSuccess);
	}

BuffOfs = sprintf(szBuff,",\"Seq\"");
for(Idx = 0; Idx < NumInputFilesProcessed; Idx++)
	BuffOfs += sprintf(&szBuff[BuffOfs],",\"%s\"",LocateFileName(Idx+1));
CUtility::SafeWrite(m_hOutFile,szBuff,BuffOfs);
BuffOfs = 0;

// construct an index of ptrs to each and initialise
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Constructing index for %d reads\n",m_NumDescrReads);

if(m_ppReadsIdx == NULL)
	{
	size_t memreq = sizeof(tsRawReadV6 *) * m_NumDescrReads;
#ifdef _WIN32
	m_ppReadsIdx = (tsRawReadV6 **) malloc((size_t)memreq);
	if(m_ppReadsIdx == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Reads Idx: Memory allocation of %lld bytes failed",(INT64)memreq);
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_ppReadsIdx = (tsRawReadV6 **)mmap(NULL,(size_t)memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_ppReadsIdx == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %lld bytes through mmap()  failed",(INT64)memreq,strerror(errno));
		m_ppReadsIdx = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_ReadsIdxBuffAllocMem = memreq;
	}

pRead = (tsRawReadV6 *)m_pDataBuff;
BuffOfs = 0;
for(Idx = 0; Idx < (int)m_NumDescrReads; Idx++)
	{
	m_ppReadsIdx[Idx] = pRead;
	BuffOfs += sizeof(tsRawReadV6) + pRead->ReadLen + pRead->DescrLen;
	pRead = (tsRawReadV6 *)&m_pDataBuff[BuffOfs];
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: sorting %d reads",m_NumDescrReads);

m_MTqsort.qsort(m_ppReadsIdx,(INT64)m_NumDescrReads,sizeof(tsRawReadV6 *),SortReads);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: sort of %d reads completed",m_NumDescrReads);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: starting to identify those reads with unique sequences");

NumUniqueReadSeqs = 0;

pRead = m_ppReadsIdx[0];
memset(SampleReadCnts,0,NumInputFilesProcessed * sizeof(UINT32));
SampleReadCnts[pRead->FileID-1] = 1;
bMinAnySampleCnts = MinSampleCnts == 1 ? true : false;
NumSamplesTotalCnts = 1;
BuffOfs = 0;
for(Idx = 1; Idx < (int)m_NumDescrReads; Idx++)
	{
	pRead1 = m_ppReadsIdx[Idx];
	Rslt = CompareRead(pRead,pRead1);
	if(Rslt > 0)							// problem if pRead ever compares as lexicographically higher than pRead1
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: Unexpected read lexicographic order");

	if(Rslt != 0)							// read is not identical to previous read
		{
		if(bMinAnySampleCnts && NumSamplesTotalCnts >= MinTotalCnts)
			{
			// output current set of counts for pRead here
			if((BuffOfs + 1000) > sizeof(szBuff))
				{
				CUtility::SafeWrite(m_hOutFile,szBuff,(unsigned int)BuffOfs);
				BuffOfs = 0;
				}

			pDst = Sequence;
			pSrc = &pRead->Read[pRead->DescrLen+1];
			for(SeqIdx =0; SeqIdx < pRead->ReadLen; SeqIdx++,pDst++,pSrc++)
				*pDst = *pSrc & 0x07;

			BuffOfs += sprintf(&szBuff[BuffOfs],"\n%d,\"%s\"",NumUniqueReadSeqs+1,CSeqTrans::MapSeq2Ascii(Sequence,pRead->ReadLen));
			for(SmplIdx = 0; SmplIdx < NumInputFilesProcessed; SmplIdx+=1)
				{
				BuffOfs += sprintf(&szBuff[BuffOfs],",%d",SampleReadCnts[SmplIdx]);
				}
			NumUniqueReadSeqs += 1;
			}
		memset(SampleReadCnts,0,NumInputFilesProcessed * sizeof(UINT32));
		pRead = pRead1;
		SampleReadCnts[pRead->FileID-1] = 1;
		NumSamplesTotalCnts = 1;
		bMinAnySampleCnts = MinSampleCnts == 1 ? true : false;
		NumUniqueReadSeqs += 1;
		continue;
		}

	// reads are identical, attribute counts to file from which that read was parsed
	SampleReadCnts[pRead1->FileID-1] += 1;
	if(!bMinAnySampleCnts && SampleReadCnts[pRead1->FileID-1] >= (UINT32)MinSampleCnts)
		bMinAnySampleCnts = true;
	NumSamplesTotalCnts += 1;

	}

if(bMinAnySampleCnts && NumSamplesTotalCnts >= MinTotalCnts)
	{
	// output current set of counts for pRead here
	pRead = pRead1;
	if((BuffOfs + 3000) > sizeof(szBuff))
		{
		CUtility::SafeWrite(m_hOutFile,szBuff,(unsigned int)BuffOfs);
		BuffOfs = 0;
		}
	pDst = Sequence;
	pSrc = &pRead->Read[pRead->DescrLen+1];
	for(SeqIdx =0; SeqIdx < pRead->ReadLen; SeqIdx++,pDst++,pSrc++)
		*pDst = *pSrc & 0x07;

	BuffOfs += sprintf(&szBuff[BuffOfs],"\n%d,\"%s\"",NumUniqueReadSeqs+1,CSeqTrans::MapSeq2Ascii(Sequence,pRead->ReadLen));
	for(SmplIdx = 0; SmplIdx < NumInputFilesProcessed; SmplIdx+=1)
		{
		BuffOfs += sprintf(&szBuff[BuffOfs],",%d",SampleReadCnts[SmplIdx]);
		}
	NumUniqueReadSeqs += 1;
	}
if(BuffOfs)
	CUtility::SafeWrite(m_hOutFile,szBuff,(unsigned int)BuffOfs);
close(m_hOutFile);
m_hOutFile = -1;

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Process: There were %d unique read sequences",NumUniqueReadSeqs);

return((teBSFrsltCodes)Rslt);
}




teBSFrsltCodes
CProcRawReads::LoadReads(bool bIsPairRead,					// true if this file to process contains the paired reads
		  UINT32 PairReadID,				// if non-zero then start paired reads identifiers from this value and increment after each read processed
		  UINT32 NumReadsLimit,				// limit to at most this number of reads, 0 if no limit
		  etFQMethod Quality,				// fastq quality value method
  		  int Trim5,						// trim this many bases off leading sequence 5' end
		  int Trim3,						// trim this many bases off trailing sequence 3' end
		  int FileID,						// uniquely identifies source file
		  char *pszFile)					// process from this file
{
static int FileNamesOfs = 0;
teBSFrsltCodes Rslt;
int Idx;
int NumDescrReads;
bool bIsFastq;
bool bSimReads;

int DescrLen;
UINT8 szDescrBuff[cRRMaxDescrLen];
int ReadLen;
UINT8 szReadBuff[cMaxReadLen];
int QualLen;
UINT8 szQualBuff[cMaxReadLen];
UINT8 *pReadBuff;
UINT8 *pQualBuff;
UINT8 Qphred;
int NumInvalValues = 0;
int NumUnsupportedBases = 0;
int NumUnderlength = 0;
unsigned long Elapsed;
unsigned long CurElapsed;
CStopWatch CurTime;
CurTime.Start();
Elapsed = CurTime.ReadUSecs();


CFasta Fasta;
if((Rslt=(teBSFrsltCodes)Fasta.Open(pszFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: Unable to open '%s' [%s] %s",pszFile,Fasta.ErrText((teBSFrsltCodes)Rslt),Fasta.GetErrMsg());
	Reset();
	return(Rslt);
	}
if(!m_FileHdr.NumFiles)
	FileNamesOfs = 0;

m_FileHdr.NumFiles += 1;
if((FileNamesOfs + strlen(pszFile) + 1) < sizeof(m_FileHdr.FileNames))
	{
	strcpy((char *)&m_FileHdr.FileNames[FileNamesOfs],pszFile);
	FileNamesOfs += (int)strlen(pszFile) + 1;
	}

m_bIsSOLiD = Fasta.IsSOLiD();

bIsFastq = Fasta.IsFastq();

NumUnsupportedBases = 0;
NumDescrReads = 0;
bSimReads = false;
Elapsed = CurTime.ReadUSecs();
while((NumReadsLimit == 0 || NumDescrReads < (int)NumReadsLimit) && ((Rslt = (teBSFrsltCodes)(ReadLen = Fasta.ReadSequence(szReadBuff,sizeof(szReadBuff),true,false))) > eBSFSuccess))
	{
	if(m_bActivity && !(m_NumDescrReads % 25000))		// reduces the overhead of checking on elapsed time if only checked every 25000 reads...
		{
		CurElapsed = CurTime.ReadUSecs();
		if((CurElapsed - Elapsed) > 60)
			{
			Elapsed = CurElapsed;
			gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing read %d",m_NumDescrReads);
			}
		}

	NumDescrReads += 1;
	if(ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		DescrLen = Fasta.ReadDescriptor((char *)szDescrBuff,sizeof(szDescrBuff));
		ReadLen = Fasta.ReadSequence(szReadBuff,sizeof(szReadBuff));

		if(ReadLen < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads processed",NumDescrReads);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szDescrBuff);
			Fasta.Close();
			Reset();
			return(eBSFerrParse);
			}

		// on first read only, check if these reads have been simulated, in which case the descriptor contains the loci of where the read was simulated from
		if(NumDescrReads == 1)
			{
			if(!strncmp((char *)szDescrBuff,"lcl|usimreads|",14))
				bSimReads = true;
			else
				bSimReads = false;
			}
		// ensure would still have a sequence of at least cRRMinSeqLen after any end trims applied
		if((Trim5 + Trim3 + (int)cRRMinSeqLen) > ReadLen)
			{
			if(!NumUnderlength++)
				{
				if((Trim5 + Trim3) > 0)
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: at least one under length ( < %d) sequence of trimmed length %d after end trims has been sloughed..",cRRMinSeqLen,ReadLen - (Trim5 + Trim3));
				else
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: at least one under length ( < %d ) sequence of length %d has been sloughed..",cRRMinSeqLen,ReadLen);
				}
			continue;
			}

		if(bIsFastq && Quality != eFQIgnore)
			{
			QualLen = Fasta.ReadQValues((char *)szQualBuff,sizeof(szQualBuff));
			if(QualLen != ReadLen)		// must be same...
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Load: quality length (%d) not same as read length (%d) for '%s' entry in file '%s'",QualLen,ReadLen,szDescrBuff,pszFile);
				Fasta.Close();
				Reset();
				return(eBSFerrParse);
				}
			// normalise the quality score to be in range 0..15 (needs to fit into 4 bits!)
			pQualBuff = szQualBuff;
			for(Idx = 0; Idx < ReadLen; Idx++,pQualBuff++)
				{
				switch(Quality) {
					case eFQIgnore:		// simply treat as the minimum phred
						Qphred = 0;
						break;

					case eFQSanger:		// Qphred = -10 log10(P), where Qphred is in range 0..93
						if(*pQualBuff < 33 || *pQualBuff >= 126)
							{
							if(!NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Sanger",*(char *)pQualBuff,szDescrBuff,pszFile);
							if(*pQualBuff < 33)
								*pQualBuff = 33;
							else
								*pQualBuff = 125;
							}
						Qphred = *pQualBuff - 33;	// Sanger encodes into ascii starting from decimal 33 '!'
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;

					case eFQIllumia:	//Qphred = -10 log10(P), where Qphred is in range 0..63
						if(*pQualBuff < 64 || *pQualBuff >= 126)
							{
							if(!NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Illumina 1.3+",*(char *)pQualBuff,szDescrBuff,pszFile);
							if(*pQualBuff < 64)
								*pQualBuff = 64;
							else
								*pQualBuff = 125;
							}

						Qphred = *pQualBuff - 64;	//Illumia encodes into ascii starting from decimal 64
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						break;

					case eFQSolexa:		// SolexaQ = -10 log10(P/(1-P)), where SolexaQ is in range -5 to 62, note the negative value
										// negative values will result if P > 0.5
										// $Q = 10 * log(1 + 10 ** (ord(SolexaQphred) - 64) / 10.0)) / log(10);
										// once Qphred is over about 15 then essentially same as Sanger and Illumina 1.3+ so
										// is it worth doing the full conversion????
						if(*pQualBuff < 59 || *pQualBuff >= 126)
							{
							if(!NumInvalValues++)
								gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: quality value of '%c' outside range expected in read '%s' from file '%s' for Solexa/Illumina pre 1.3",*(char *)pQualBuff,szDescrBuff,pszFile);
							if(*pQualBuff < 64)
								*pQualBuff = 64;
							else
								*pQualBuff = 125;
							}
						Qphred = *pQualBuff - 59;	// Solexa/Illumina encodes into ascii starting from decimal 59
						Qphred = (UINT8)(10 * log(1 + pow(10.0,((double)Qphred/10.0) / log(10.0))));	//
						if(Qphred > 40)				// clamp at phred equiv to 0.0001
							Qphred = 40;
						*pQualBuff = (UINT8)(((UINT32)Qphred * 31) / 40);
						break;
					}
				*pQualBuff = (UINT8)((((UINT32)Qphred+2)*15)/40);
				}
			// pack the read and quality, read into the low order bits 0..3, quality into bits 4..7
			pQualBuff = szQualBuff;
			pReadBuff = szReadBuff;
			for(Idx = 0; Idx < ReadLen; Idx++,pQualBuff++,pReadBuff++)
				szReadBuff[Idx] |= *pQualBuff << 4;
			}

		// apply any end trims
		if(Trim5)
			{
			ReadLen -= Trim5;
			memmove(szReadBuff,&szReadBuff[Trim5],ReadLen);
			}
		if(Trim3)
			ReadLen -= Trim3;

		// truncate descriptors at 1st whitespace unless fasta was generated by simulating reads in which case
		// the whole descriptor is retained as where the read was simulated from is of interest
		if(!bSimReads)	// if not simulated reads
			{
			for(Idx = 0; Idx < cMaxDescrIDLen-1; Idx++)
				{
				if(szDescrBuff[Idx] == '\0' || isspace(szDescrBuff[Idx]))
					break;
				}
			szDescrBuff[Idx] = '\0';
			DescrLen = Idx;
			}
		if((Rslt=(teBSFrsltCodes)AddEntry(bIsPairRead,PairReadID,FileID,DescrLen,(char *)szDescrBuff,ReadLen,szReadBuff))!=eBSFSuccess)
			break;
		if(PairReadID > 0)
			PairReadID += 1;
		continue;
		}
	else
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszFile,
				Rslt > 0 ? "overlength read" : "not a multifasta short reads or fastq file");
		Fasta.Close();
		Reset();
		return(eBSFerrParse);
		}
	}
if(Rslt != eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Errors processing file: %s ",pszFile);
	while(Fasta.NumErrMsgs())
		gDiagnostics.DiagOut(eDLFatal,gszProcName,Fasta.GetErrMsg());
	Fasta.Close();
	Reset();
	return(Rslt);
	}
if(m_bActivity)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Processing read %d",m_NumDescrReads);
Fasta.Close();
if(NumInvalValues > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unexpected quality values read from file '%s'",NumInvalValues,pszFile);
if(NumUnsupportedBases > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d unsupported bases read from file '%s'",NumUnsupportedBases,pszFile);
if(NumUnderlength > 0)
	gDiagnostics.DiagOut(eDLWarn,gszProcName,"Load: total of %d under length sequences sloughed from file '%s'",NumUnderlength,pszFile);

return(eBSFSuccess);
}



int
CProcRawReads::AddEntry(bool bIsPairRead,		// true if this is the paired read
		 UINT32 PairReadID,		// identifies partner of this read if paired read processing
		 UINT8 FileID,			// identifies file from which this read was parsed
		 int DescrLen,			// length of following descriptor
		 char *pszReadDescr,	// copy of descriptor, used to pair reads with matching descriptors
		 int ReadLen,			// length of following read
		 UINT8 *pszReadBuff)	// packed read + phred score
{
UINT8 *pTmpAlloc;
tsRawReadV6 DescrRead;
size_t memreq;
if(m_pDataBuff == NULL)
	{
	memreq = cRRDataBuffAlloc * 2;
#ifdef _WIN32
	m_pDataBuff = (UINT8 *) malloc((size_t)memreq);
	if(m_pDataBuff == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %lld bytes failed",(INT64)memreq);
		Reset();
		return(eBSFerrMem);
		}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pDataBuff = (UINT8 *)mmap(NULL,(size_t)memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pDataBuff == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory allocation of %lld bytes through mmap()  failed",(INT64)memreq,strerror(errno));
		m_pDataBuff = NULL;
		Reset();
		return(eBSFerrMem);
		}
#endif
	m_DataBuffAllocMem = memreq;
	m_DataBuffOfs = 0;
	}

// is there a need to allocate more memory? NOTE: allow safety margin of 100K bytes
if((m_DataBuffAllocMem - m_DataBuffOfs) < (sizeof(tsRawReadV6) +  ReadLen + DescrLen + 100000))
	{
	memreq = m_DataBuffAllocMem + cRRDataBuffAlloc;
#ifdef _WIN32
	pTmpAlloc = (UINT8 *) realloc(m_pDataBuff,memreq);
#else
	pTmpAlloc = (UINT8 *)mremap(m_pDataBuff,m_DataBuffAllocMem,memreq,MREMAP_MAYMOVE);
	if(pTmpAlloc == MAP_FAILED)
		pTmpAlloc = NULL;
#endif
	if(pTmpAlloc == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddEntry: Memory reallocation to %lld bytes failed - %s",memreq,strerror(errno));
		Reset();
		return(eBSFerrMem);
		}
	m_pDataBuff = pTmpAlloc;
	m_DataBuffAllocMem = memreq;
	}
DescrRead.ReadID = ++m_NumDescrReads;
DescrRead.PairReadID = PairReadID;
if(bIsPairRead)
	DescrRead.PairReadID |= 0x80000000;
DescrRead.ReadLen = ReadLen;
DescrRead.NumReads = 1;
DescrRead.DescrLen=(UINT8)DescrLen;
DescrRead.FileID = FileID;
memmove(&m_pDataBuff[m_DataBuffOfs],&DescrRead,sizeof(DescrRead)-1);
m_DataBuffOfs += (int)sizeof(DescrRead)-1;
memmove(&m_pDataBuff[m_DataBuffOfs],pszReadDescr,DescrLen+1);
m_DataBuffOfs += DescrRead.DescrLen + 1;
memmove(&m_pDataBuff[m_DataBuffOfs],pszReadBuff,ReadLen);
m_DataBuffOfs += ReadLen;
return(eBSFSuccess);
}


// CompareRead
// Returns 0 if reads are identical in both sequence and length
// Returns -1 if pRead1 is lexigraphically ordered less than pRead2 or if same then shorter length
// Returns +1 if pRead2 is lexographically ordered higher than pRead2 or if same then longer length
int
CProcRawReads::CompareRead(tsRawReadV6 *pRead1,tsRawReadV6 *pRead2)
{
int Idx;
UINT8 *pSeq1 = &pRead1->Read[pRead1->DescrLen+1];
UINT8 *pSeq2 = &pRead2->Read[pRead2->DescrLen+1];

for(Idx = 0; Idx < pRead1->ReadLen; Idx++, pSeq1++, pSeq2++)
	{
	if(Idx >= pRead2->ReadLen)
		return(1);
	if((*pSeq1 & 0x07) < (*pSeq2 & 0x07))
		return(-1);
	if((*pSeq1 & 0x07) > (*pSeq2 & 0x07))
		return(1);
	}
if(pRead1->ReadLen < pRead2->ReadLen)
	return(-1);
return(0);
}

// SortReads
// Used by qsort() for indexing read sequences into suffix array
int
CProcRawReads::SortReads(const void *arg1, const void *arg2)
{
int Idx;
UINT8 *pSeq1;
UINT8 *pSeq2;

tsRawReadV6 *pEl1 = *(tsRawReadV6 **)arg1;
tsRawReadV6 *pEl2 = *(tsRawReadV6 **)arg2;

pSeq1 = &pEl1->Read[pEl1->DescrLen+1];
pSeq2 = &pEl2->Read[pEl2->DescrLen+1];

for(Idx = 0; Idx < pEl1->ReadLen; Idx++, pSeq1++, pSeq2++)
	{
	if(Idx >= pEl2->ReadLen)
		return(1);
	if((*pSeq1 & 0x07) < (*pSeq2 & 0x07))
		return(-1);
	if((*pSeq1 & 0x07) > (*pSeq2 & 0x07))
		return(1);
	}
if(pEl1->ReadLen < pEl2->ReadLen)
	return(-1);
return(0);
}

// SortReadsPair
// Sort reads pairs by descriptor and file ascending
int
CProcRawReads::SortReadsPair(const void *arg1, const void *arg2)
{
int Rslt;

tsRawReadV6 *pEl1 = *(tsRawReadV6 **)arg1;
tsRawReadV6 *pEl2 = *(tsRawReadV6 **)arg2;

if((Rslt = stricmp((char *)pEl1->Read,(char *)pEl2->Read))!=0)
	return(Rslt);
return(pEl1->FileID - pEl2->FileID);
}
