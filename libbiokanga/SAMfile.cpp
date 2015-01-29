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

// BAM file format
// http://samtools.sourceforge.net/SAMv1.pdf
// https://github.com/samtools/hts-specs
// ---------------
// Header:
// magic	| char[4]		| BAM magic string BAM\1
// l_text	| int32_t		| text Length of the header text, including any NULL padding
// text		| char[l_text]	| Plain header text in SAM; not necessarily NULL terminated | char[l text]
// n_refs	| int32_t		| ref # reference sequences 
//
// Repeat following reference sequences for n_refs:
// l_name	| int32_t		| name Length of the reference name plus 1 (including NULL)
// name		| char[l_name]	| Reference sequence name; NULL terminated
// l ref	| int32_t       | Length of the reference sequence
//     
// Repeat following alignments:
// block_size | int32        | Length of the remainder of this alignment record (includes any auxiliary data)
// refID	  | int32        | Reference sequence ID,  -1 <= refID < n_ref; -1 for a read without a mapping position
// pos		  | int32        | 0-based leftmost coordinate (= POS - 1) 
// bin_mq_nl  | uint32       | bin<<16|MAPQ<<8|l_read_name ; bin is computed by the reg2bin(); l_read_name is the length of read name below (= length(QNAME) + 1).
// flag_nc    | uint32       | FLAG<<16|n_cigar_op; n_cigar_op is the number of operations in CIGAR
// l_seq      | int32        | Length of SEQ
// next_refID | int32        | RefID of the next segment (-1 <= mate_refID < n_ref)
// next_pos   | int32        |  0-based leftmost pos of the next segment (= PNEXT - 1)
// tlen       | int32        | Template length (= TLEN)
// read_name  | char[l_read name] | NULL terminated (QNAME plus a tailing `\0')
// cigar      | uint32[n_cigar_op] | CIGAR: op_len<<4|op. `MIDNSHP=X'!`012345678' 
// seq        | uint8 t[(l_seq+1)/2]| 4-bit encoded read: `=ACMGRSVTWYHKDBN'! [0; 15]; other characters mapped to `N'; high nybble firrst (1st base in the highest 4-bit of the 1st byte)
// qual       | char[l_seq]  | Phred base quality (a sequence of 0xFF if absent)
//
//     Repeat auxiliary data as may be required:
// tag        | char[2]      | Two-character tag_char
// val_type   | char         | Value type: AcCsSiIfZHB
// value	  | ?????        | Tag value (by val type)

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include <sys/mman.h>
#include "./commhdrs.h"
#endif
#include "./bgzf.h"

const char m_CigarOpsMap[] = {'M','I','D','N','S','H','P','=','X'};
const char m_BasesMap[] = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};

CSAMfile::CSAMfile(void)
{
m_pBGZF = NULL;
m_hOutSAMfile = -1;
m_hOutBAIfile = -1;
m_gzOutSAMfile = NULL;
m_pgzOutCSIfile = NULL;
m_hInSAMfile = -1;
m_gzInSAMfile = NULL;
m_pInBGZF = NULL;
m_pBAM = NULL;
m_pBAI = NULL;
m_pRefSeqs = NULL;
m_pBAIChunks = NULL;
m_pChunkBins = NULL;
m_p16KOfsVirtAddrs = NULL;
Reset(false);
}


CSAMfile::~CSAMfile(void)
{
Reset(false);
}

void
CSAMfile::Reset(bool bSync) // if bSync true then fsync before closing output file handles
{

if(m_gzOutSAMfile != NULL)
	{
	gzclose(m_gzOutSAMfile);
	m_gzOutSAMfile = NULL;
	}

if(m_gzInSAMfile != NULL)
	{
	gzclose(m_gzInSAMfile);
	m_gzInSAMfile = NULL;
	}

if(m_hOutSAMfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutSAMfile);
#else
		fsync(m_hOutSAMfile);
#endif
	close(m_hOutSAMfile);
	m_hOutSAMfile = -1;
	}

if(m_hInSAMfile != -1)
	{
	close(m_hInSAMfile);
	m_hInSAMfile = -1;
	}

if(m_hOutBAIfile != -1)
	{
	if(bSync)
#ifdef _WIN32
		_commit(m_hOutBAIfile);
#else
		fsync(m_hOutBAIfile);
#endif
	close(m_hOutBAIfile);
	m_hOutBAIfile = -1;
	}

if(m_pgzOutCSIfile != NULL)
	{
	if(bSync)
		bgzf_flush(m_pgzOutCSIfile);
	bgzf_close(m_pgzOutCSIfile);
	m_pgzOutCSIfile = NULL;
	}

if(m_pBGZF != NULL)
	{
	if(bSync)
		bgzf_flush(m_pBGZF);
	bgzf_close(m_pBGZF);
	m_pBGZF = NULL;
	}

if(m_pInBGZF != NULL)
	{
	bgzf_close(m_pInBGZF);
	m_pInBGZF = NULL;
	}

if(m_pBAM != NULL)
	{
#ifdef _WIN32
	free(m_pBAM);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pBAM != MAP_FAILED)
		munmap(m_pBAM,m_AllocBAMSize);
#endif
	m_pBAM = NULL;
	}

if(m_pBAI != NULL)
	{
#ifdef _WIN32
	free(m_pBAI);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pBAI != MAP_FAILED)
		munmap(m_pBAI,m_AllocBAISize);
#endif
	m_pBAI = NULL;
	}

if(m_pBAIChunks != NULL)
	{
#ifdef _WIN32
	free(m_pBAIChunks);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pBAIChunks != MAP_FAILED)
		munmap(m_pBAIChunks,m_AllocBAIChunks * sizeof(tsBAIChunk));
#endif
	m_pBAIChunks = NULL;
	}


if(m_pChunkBins != NULL)
	{
#ifdef _WIN32
	free(m_pChunkBins);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pChunkBins != MAP_FAILED)
		munmap(m_pChunkBins,cNumSAIBins * sizeof(tsBAIbin));
#endif
	m_pChunkBins = NULL;
	}

if(m_p16KOfsVirtAddrs != NULL)
	{
#ifdef _WIN32
	free(m_p16KOfsVirtAddrs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_p16KOfsVirtAddrs != MAP_FAILED)
		munmap(m_p16KOfsVirtAddrs,m_Alloc16KOfsVirtAddrsSize);
#endif
	m_p16KOfsVirtAddrs = NULL;
	}

if(m_pRefSeqs != NULL)
	{
#ifdef _WIN32
	free(m_pRefSeqs);				// was allocated with malloc/realloc, or mmap/mremap, not c++'s new....
#else
	if(m_pRefSeqs != MAP_FAILED)
		munmap(m_pRefSeqs,m_AllocRefSeqsSize);
#endif
	m_pRefSeqs = NULL;
	}

m_ComprLev = 0;
m_AllocBAMSize = 0;
m_CurBAMLen = 0;
m_NumBAMSeqNames = 0;
m_CurInBAMIdx = 0;
m_bInEOF = false;
m_TotInBAMProc = 0;
m_InBAMHdrLen = 0;

m_AllocBAISize = 0;
m_CurBAILen = 0;
m_NumBAISeqNames = 0;

m_NumBinsWithChunks = 0;
m_AllocBAIChunks = 0;
m_NumChunks = 0;
m_NumOf16Kbps = 0;
m_Alloc16KOfsVirtAddrsSize = 0;
m_MaxIdxRefSeqLen = 0;
m_MaxRefSeqLen = 0;

m_CSI_min_shift = 14;
m_CSI_depth = 5;

m_AllocRefSeqsSize = 0;
m_CurRefSeqsLen = 0;
m_NumRefSeqNames = 0;
m_CurRefSeqNameID = 0;

m_bBAMfile = false;	
m_szSAMfileName[0] = '\0';
m_szBAIfileName[0] = '\0';

m_ParseSeqState = 0;
m_szParsedDescriptor[0] = '\0';
m_ParsedDescrLen = 0;
m_szParsedSeqBases[0] = '\0';
m_ParsedSeqLen = 0;

m_LocateRefSeqHistDepth = 0;
m_szLastNotLocatedRefSeqName[0] = '\0';
}

bool									// open and check if a SAM or BAM format file, returns true if was a SAM or BAM
CSAMfile::IsSAM(char *pszSAMFile)	// expected to be a SAM(gz) or if extension '.BAM' then a BAM file
{
int FileNameLen;
int hInSAMfile;
eSAMFileType SAMType;
gzFile gzInSAMfile;
char szBAMHdr[1000];
BGZF* pInBGZF;						// BAM is BGZF compressed
int CurBAMLen;


if(pszSAMFile == NULL || pszSAMFile[0] == '\0')
	return(false);
// classify file type by extension
// if file has suffix of '.bam' then classify as eSFTBAM
// if file has suffix of '.gz' then classify as eSFTSAMgz
// if file has any other extension then classify as SAM
// 
FileNameLen = (int)strlen(pszSAMFile);
if(FileNameLen > 3 && !stricmp(&pszSAMFile[FileNameLen-3],".gz"))
   SAMType = eSFTSAMgz;
else
	{
	if(FileNameLen > 4 && !stricmp(&pszSAMFile[FileNameLen-4],".bam"))
		SAMType = eSFTBAM;
	else
		SAMType = eSFTSAM;
	}
if(SAMType != eSFTSAMgz)
	{
#ifdef _WIN32
	hInSAMfile = open(pszSAMFile,( O_RDONLY | _O_BINARY | _O_SEQUENTIAL),_S_IREAD);
#else
	hInSAMfile = open(pszSAMFile,O_RDONLY,S_IREAD);
#endif
	if(hInSAMfile < 0)
		return(false);			
	}
else
	{
	gzInSAMfile = gzopen(pszSAMFile,"rb");
	if(gzInSAMfile == NULL)
		return(false);	
	}

if(SAMType >= eSFTBAM)
	{
	// BAM will using BGZF compression ..
	if((pInBGZF = bgzf_dopen(hInSAMfile, "r"))==NULL)
		return(false);	
	hInSAMfile = -1;

	// try reading the header, bgzf_read will confirm it does start with "BAM\1" ....
	if((CurBAMLen = (int)bgzf_read(pInBGZF,szBAMHdr,100)) < 100)		// will be -1 if errors ...
		{
		bgzf_close(pInBGZF);
		return(false);	
		}
	bgzf_close(pInBGZF);
	}
else
	{
	// try reading in intial header and check that it does look like a SAM
	// accepting as SAM if first line starts with '@HD\tVN:'
	if(SAMType == eSFTSAMgz)
		{
		if((CurBAMLen = gzread(gzInSAMfile,szBAMHdr,100)) < 100)
			{
			gzclose(gzInSAMfile);
			return(false);
			}
		gzclose(gzInSAMfile);
		}
	else
		{
		if((CurBAMLen = read(hInSAMfile,szBAMHdr,100)) < 100)
			{
			close(hInSAMfile);
			return(false);
			}
		close(hInSAMfile);
		}
	if(strnicmp((const char *)szBAMHdr,"@HD\tVN:",6))
		return(false);
	}

return(true);
}



int					// open and initiate processing for SAM/BAM reads processing
CSAMfile::Open(char *pszSAMFile)	// SAM(gz) or BAM file name
{
int FileNameLen;
eSAMFileType SAMType;
if(pszSAMFile == NULL || pszSAMFile[0] == '\0')
	return(eBSFerrParams);

Reset();

// classify file type by extension
// if file has suffix of '.bam' then classify as eSFTBAM
// if file has suffix of '.gz' then classify as eSFTSAMgz
// if file has any other extension then classify as SAM
// 
FileNameLen = (int)strlen(pszSAMFile);
if(FileNameLen > 3 && !stricmp(&pszSAMFile[FileNameLen-3],".gz"))
   SAMType = eSFTSAMgz;
else
	{
	if(FileNameLen > 4 && !stricmp(&pszSAMFile[FileNameLen-4],".bam"))
		SAMType = eSFTBAM;
	else
		SAMType = eSFTSAM;
	}

m_SAMFileType = SAMType;
strcpy(m_szSAMfileName,pszSAMFile);

if(SAMType != eSFTSAMgz)
	{
#ifdef _WIN32
	m_hInSAMfile = open(m_szSAMfileName,( O_RDONLY | _O_BINARY | _O_SEQUENTIAL),_S_IREAD);
#else
	m_hInSAMfile = open(m_szSAMfileName,O_RDONLY,S_IREAD);
#endif
	if(m_hInSAMfile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: unable to open for reading file '%s'",m_szSAMfileName);
		Reset();
		return(eBSFerrOpnFile);
		}
	}
else
	{
	m_gzInSAMfile = gzopen(m_szSAMfileName,"rb");
	if(m_gzInSAMfile == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: unable to open for reading gzip'd file '%s'",m_szSAMfileName);
		return(eBSFerrOpnFile);
		}
	gzbuffer(m_gzInSAMfile,cAllocSAMSize);		// buffering to reduce number of disk reads required
	}

// make initial mem allocation
if(SAMType == eSFTSAM || SAMType == eSFTSAMgz)
	m_AllocBAMSize = cAllocSAMSize;
else
	m_AllocBAMSize = cAllocBAMSize;

if((m_pBAM = (UINT8 *)malloc(m_AllocBAMSize))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to alloc memory for bufferingr");
	Reset();
	return(eBSFerrMem);
	}

if((m_pRefSeqs = (tsRefSeq *)malloc(cAllocRefSeqSize))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to alloc memory for reference sequence names");
	Reset();
	return(eBSFerrMem);
	}
m_AllocRefSeqsSize = cAllocRefSeqSize;
m_NumRefSeqNames = 0;
m_CurRefSeqsLen = 0;

// for BAM then alloc to hold reference sequence names
if(SAMType >= eSFTBAM)
	{
	// BAM will using BGZF compression ..
	if((m_pInBGZF = bgzf_dopen(m_hInSAMfile, "r"))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: unable to initialise for BGZF processing on file '%s'",m_szSAMfileName);
		Reset();
		return(eBSFerrOpnFile);
		}
	m_hInSAMfile = -1;

	// try reading the header, bgzf_read will confirm it does start with "BAM\1" ....
	if((m_CurBAMLen = (int)bgzf_read(m_pInBGZF,m_pBAM,100)) < 100)		// will be -1 if errors ...
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: Not a BAM format file '%s'",m_szSAMfileName);
		Reset();
		return(eBSFerrOpnFile);
		}
	m_bInEOF = false;
	m_CurBAMLen = 100;
	m_CurInBAMIdx = 4;
	m_TotInBAMProc = 4;
	m_InBAMHdrLen = 0;
	}
else
	{
	// try reading in intial header and check that it does look like a SAM
	// accepting as SAM if first line starts with '@HD\tVN:'
	if(SAMType == eSFTSAMgz)
		{
		if((m_CurBAMLen = gzread(m_gzInSAMfile,m_pBAM,100)) < 100)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: Not a SAM format file '%s'",m_szSAMfileName);
			Reset();
			return(eBSFerrOpnFile);
			}
		gzseek(m_gzInSAMfile,0,SEEK_SET);
		}
	else
		{
		if((m_CurBAMLen = read(m_hInSAMfile,m_pBAM,100)) < 100)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: Not a SAM format file '%s'",m_szSAMfileName);
			Reset();
			return(eBSFerrOpnFile);
			}
		lseek(m_hInSAMfile,0,SEEK_SET);
		}
	if(strnicmp((const char *)m_pBAM,"@HD\tVN:",6))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Open: Not a SAM format file '%s'",m_szSAMfileName);
		Reset();
		return(eBSFerrOpnFile);
		}
	m_bInEOF = false;
	m_CurInBAMIdx = 0;
	m_TotInBAMProc = 0;
	}

// initialisation completed
return(eBSFSuccess);
}


char *
CSAMfile::TrimWhitespace(char *pTxt)
{
char *pStart;
char Chr;
	// strip leading whitespace
while(Chr = *pTxt++)
	if(!isspace(Chr))
			break;
if(Chr == '\0')					// empty line?
	return(pTxt-1);
pStart = pTxt-1;
while(Chr = *pTxt)			// fast forward to line terminator
	pTxt++;
pTxt-=1;
while(Chr = *pTxt--)
	if(!isspace(Chr))
		break;
pTxt[2] = '\0';
return(pStart);
}

UINT32										// returns estimated number of sequences in SAM or BAM file
CSAMfile::EstSizes(char *pszFile,			// SAM or BAM file path+name to estimate sizes
			  INT64 *pFileSize,				// file is this size on disk
			  INT32 *pEstMaxDescrLen,		// with estimated maximum descriptor length
			  INT32 *pEstMeanDescrLen,		// estimated mean descriptor length
			  INT32 *pEstMaxSeqLen,			// and estimated maximum sequence length
			  INT32 *pEstMeanSeqLen,		// estimated mean sequence length
			  INT32 *pEstScoreSchema)		// currently will always return 0: no scoring 
{
int Rslt;
int NumSeqs;
UINT32 EstTotSeqs;
INT64 SumDescrLens;
INT64 SumSeqLens;
int MaxDescrLen;
int MinDescrLen;
int MaxSeqLen;
int MinSeqLen;
int MeanDescrLen;
int MeanSeqLen;
INT64 FileSize;

NumSeqs = 0;
SumDescrLens = 0;
SumSeqLens = 0;
MaxDescrLen = 0;
MinDescrLen = 0;
MaxSeqLen = 0;
MinSeqLen = 0;
MeanDescrLen = 0;
MeanSeqLen = 0;
EstTotSeqs = 0;

if(pEstMaxDescrLen != NULL)
	*pEstMaxDescrLen = 0;
if(pEstMeanDescrLen != NULL)
	*pEstMeanDescrLen = 0;
if(pEstMaxSeqLen != NULL)
	*pEstMaxSeqLen = 0;
if(pEstMeanSeqLen != NULL)
	*pEstMeanSeqLen = 0;
if(pEstScoreSchema != NULL)
	*pEstScoreSchema = 0;


#ifdef _WIN32
struct _stat64 st;
if(!_stat64(pszFile,&st))
#else
struct stat64 st;
if(!stat64(pszFile,&st))
#endif
	FileSize = (INT64)st.st_size;
else
	FileSize = 0;

if(pFileSize != NULL)
	*pFileSize = FileSize;

if(FileSize < 100 || (Rslt = Open(pszFile)) != eBSFSuccess)
	return(0);

// scale up FileSize according to guestimated compression ratios
switch(m_SAMFileType) {
	case eSFTSAM:			// SAM raw text file, no compression but still have header overhead
		break;

	case eSFTSAMgz:			// SAM which was compressed with gzip
		FileSize *= 3;  
		break;

	case eSFTBAM:			// BAM bgzf compressed file
	case eSFTBAM_BAI:		// BAM bgzf plus associated BAI file
	case eSFTBAM_CSI:		// BAM bgzf plus associated CSI file
		FileSize *= 4;
		break;
	}

// read upto 100K alignments to get estimates of descriptor and sequence lengths
for(NumSeqs = 0; NumSeqs < 100000; NumSeqs++)
	{
	if((ReadSequence(NULL,0,false,false)) != eBSFFastaDescr)
		break;

	if(MaxDescrLen == 0)
		MinDescrLen = m_ParsedDescrLen;
	if(m_ParsedDescrLen > MaxDescrLen)
		MaxDescrLen = m_ParsedDescrLen;
	if(m_ParsedDescrLen < MinDescrLen)
		MinDescrLen = m_ParsedDescrLen;

	if(MaxSeqLen == 0)
		MinSeqLen = m_ParsedSeqLen;
	if(m_ParsedSeqLen > MaxSeqLen)
		MaxSeqLen = m_ParsedSeqLen;
	if(m_ParsedSeqLen < MinSeqLen)
		MinSeqLen = m_ParsedSeqLen;

	SumDescrLens += m_ParsedDescrLen;
	SumSeqLens += m_ParsedSeqLen;
	m_ParseSeqState = 0;
	}
Close();
if(NumSeqs < 1)
	return(0);

MeanDescrLen = (int)((SumDescrLens+NumSeqs-1) / NumSeqs);
MeanSeqLen = (int)((SumSeqLens+NumSeqs-1) / NumSeqs);

if(pEstMaxDescrLen != NULL)
	*pEstMaxDescrLen = MaxDescrLen;
if(pEstMeanDescrLen != NULL)
	*pEstMeanDescrLen = MeanDescrLen;
if(pEstMaxSeqLen != NULL)
	*pEstMaxSeqLen = MaxSeqLen;
if(pEstMeanSeqLen != NULL)
	*pEstMeanSeqLen = MeanSeqLen;

// estimate of total number of sequences is purely a guestimate as no idea of the compression ratios, or header section size!!!
// guestimate is that the header section will be 2% of the filesize, and that descriptor + sequences are about 75% of each alignment line

if(NumSeqs < 100000)
	EstTotSeqs = (UINT32)NumSeqs;
else
	EstTotSeqs = (UINT32)((((UINT64)FileSize * 98) / 100)  / (UINT64)(((MeanDescrLen + MeanSeqLen) * 3) / 2));

return(EstTotSeqs);
}

int 
CSAMfile::ReadDescriptor(char *pszDescriptor,int MaxLen) // copies last descriptor processed into pszDescriptor and returns copied length
{
int RetLen;
if(pszDescriptor == NULL)
	return(0);
if(m_szParsedDescriptor[0] == '\0' || m_ParsedDescrLen == 0)
	{
	pszDescriptor[0] = '\0';
	return(0);
	}
RetLen = min(MaxLen,m_ParsedDescrLen+1);
strncpy(pszDescriptor,m_szParsedDescriptor,RetLen);
pszDescriptor[RetLen-1] = '\0';
return(RetLen-1);
}

// returns next sequence in SAM/BAM file
int						// returns actual number read (eBSFSuccess == EOF,eBSFFastaDescr == End of current sequence, descriptor line now available)
CSAMfile::ReadSequence(void *pRetSeq,	// where to return sequence, can be NULL if only interested in the sequence length
					 int Max2Ret,		// max to return, ignored if pRetSeq is NULL
					 bool bSeqBase,		// if false then return ascii, if true then return as etSeqBase
					 bool RptMskUpperCase)	// default is false, UCSC softmasked use lowercase when repeat masked
{
int LineLen;
char *pTxt;
char szLine[cMaxReadLen  * 3];				// buffer input lines
etSeqBase SeqBases[cMaxReadLen*3];
char szCigar[128];
char szRNext[128];
int MAPQ;
int PNext;
int TLen;
int SeqRetLen;

switch(m_ParseSeqState) {
	case 1:					// expecting the sequence to be returned
		if(pRetSeq == NULL)
			{
			m_ParseSeqState = 0;		// next alignment is to loaded on next call to this function
			return(m_ParsedSeqLen);
			}
		SeqRetLen = min(m_ParsedSeqLen,Max2Ret);
		if(bSeqBase)
			{
			CSeqTrans::MapAscii2Sense(m_szParsedSeqBases,0,SeqBases);
			SeqBases[SeqRetLen] = eBaseEOS;
			memcpy(pRetSeq,SeqBases,SeqRetLen);
			}
		else
			{
			m_szParsedSeqBases[SeqRetLen] = '\0';
			strcpy((char *)pRetSeq,m_szParsedSeqBases);
			}
		m_ParseSeqState = 0;		// next alignment is to loaded on next call to this function
		return(SeqRetLen);
		break;

	default:
		break;
	}

m_ParseSeqState = 0;
m_szParsedDescriptor[0] = '\0';
m_ParsedDescrLen = 0;
m_szParsedSeqBases[0] = '\0';
m_ParsedSeqLen = 0;
m_ParsedFlags = 0;
m_ParsedszChrom[0] = '\0';
m_ParsedStartLoci = 0;

while((LineLen = GetNxtSAMline(szLine)) > 0)
	{
	szLine[sizeof(szLine)-1] = '\0';
	pTxt = TrimWhitespace(szLine);
	if(*pTxt=='\0' || *pTxt=='@')	// simply slough lines which were just whitespace or start with '@'
		continue;					// only interested in the alignments
	
	sscanf(szLine,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%16383s\t",m_szParsedDescriptor, &m_ParsedFlags, m_ParsedszChrom, &m_ParsedStartLoci,&MAPQ,szCigar,szRNext,&PNext,&TLen,m_szParsedSeqBases);
	m_ParsedDescrLen = (int)strlen(m_szParsedDescriptor);
	m_ParsedSeqLen = (int)strlen(m_szParsedSeqBases);
	m_ParseSeqState = 1;
	return(eBSFFastaDescr);				// descriptor is ready to be returned
	}
return(LineLen);
}


// It is expected that sequence names will have been added with AddSeqName() in alpha ascending order such that when locating
// names then the exhustative search can early terminate
int				// locates reference sequence name and returns it's SeqID, returns 0 if unable to locate a match
CSAMfile::LocateRefSeqID(char *pszRefSeqName) // reference sequence name to locate
{
int CmpRslt;
int Idx;
int NameLen;

tsRefSeq *pRefSeq; 
if(pszRefSeqName == NULL || pszRefSeqName[0] == '\0' || m_pRefSeqs == NULL || m_NumBAMSeqNames < 1)
	return(0);

// linear search; shouldn't take too long as not expecting too many reference chroms
// Wheat is a little different!
NameLen = (int)strlen(pszRefSeqName);

// one optimisation is that a short history is maintained containing the last 10 successful searches and this history is
// searched before the full reference sequence names
if(m_LocateRefSeqHistDepth)
	{
	for(Idx = 0; Idx < (int)m_LocateRefSeqHistDepth; Idx++)
		{
		pRefSeq = m_pLocateRefSeqHist[Idx];
		if(NameLen == pRefSeq->SeqNameLen && !stricmp(pszRefSeqName,pRefSeq->szSeqName))
			{
			if(Idx > 0)
				{
				for(int Idy = Idx; Idy > 0; Idy--)
					m_pLocateRefSeqHist[Idy] = m_pLocateRefSeqHist[Idy-1];
				m_pLocateRefSeqHist[0] = pRefSeq;
				}
			return(pRefSeq->SeqID);
			}
		}
	}

// no match in history so need to do a linear search
pRefSeq = m_pRefSeqs;
for(Idx = 0; Idx < (int)m_NumBAMSeqNames; Idx++)
	{
	CmpRslt = 0;
	if(NameLen == pRefSeq->SeqNameLen && (CmpRslt = stricmp(pszRefSeqName,pRefSeq->szSeqName)) == 0)
		{
		if(m_LocateRefSeqHistDepth < cMaxLocateRefSeqHist)
			m_LocateRefSeqHistDepth += 1;
		if(m_LocateRefSeqHistDepth > 1)
			{	
			for(Idx = m_LocateRefSeqHistDepth - 1; Idx > 0; Idx--)
				m_pLocateRefSeqHist[Idx] = m_pLocateRefSeqHist[Idx-1];
			}
		m_pLocateRefSeqHist[0] = pRefSeq;
		return(pRefSeq->SeqID);
		}
	pRefSeq = (tsRefSeq *)((UINT8 *)pRefSeq + pRefSeq->SeqNameLen + sizeof(tsRefSeq));
	}
return(0);
}


int				// alignment length as calculated from SAM/BAM CIGAR string, only 'M','X','=' lengths contribute
CSAMfile::CigarAlignLen(char *pszCigar)	// alignment length as calculated from SAM/BAM CIGAR
{
char *pCigar;
char CigarChr;
int CigarOPlen;
int CigarIdx;
int AlignLen;

if(pszCigar == NULL || pszCigar[0] == '\0' || pszCigar[0] == '*')
	return(0);

CigarIdx = 0;
pCigar = pszCigar;
CigarOPlen = 0;
AlignLen = 0;
while((CigarChr = *pCigar++) != '\0')
	{
	if(CigarIdx >= 50)	// can only handle CIGARs of up to 40 chars
		return(eBSFerrParse);
	if(CigarChr >= '0' && CigarChr <= '9')
		{
		CigarOPlen *= 10;
		CigarOPlen += (int)(CigarChr - '0');
		continue;
		}
	switch(tolower(CigarChr)) {
		case 'm':				// number of bases aligned (could include non-matching)
		case '=':               // number of bases aligned, exactly matching with no substitutions
		case 'x':				// number of bases aligned, none exactly matching, all substitutions
			AlignLen += CigarOPlen;
			break;

		case 'n':               // skipped region from the reference
		case 's':				// softclipped
		case 'h':				// hardclipped
		case 'p':	            // padding (silent deletion from padded reference)			
		default:				// unsupported!
			break;
		}
	if(CigarIdx > 50)
		return(eBSFerrParse);
	CigarOPlen = 0;
	}
return(AlignLen);	
}

// parses SAM formated line into a tsBAMalign structure
int											// negative if errors parsing otherwise 0 for success
CSAMfile::ParseSAM2BAMalign(char *pszSAMline, // parsing this SAM format line
				tsBAMalign *pBAMalign,        // into this tsBAMalign structure
				CBEDfile *pBEDremapper)		  // with optional remapping from features (contigs) in this BED file
{
char szDescriptor[128];			// parsed out descriptor
int Flags;						// parsed out flags
char szChrom[128];				// parsed out chrom
int ContigID;
int RelStartLoci;
int RelMateStartLoci;
int StartLoci;					
char szCigar[128];
char szRNext[128];
int MAPQ;
int PNext;
int TLen;
int LineRemainIdx;
int NumEls;
char *pszSeq;

memset(pBAMalign,0,sizeof(tsBAMalign));
if(pszSAMline == NULL || pszSAMline[0] == '\0')
	return(eBSFerrParams);

// expecting to parse as "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t", szDescriptor, Flags, m_szSAMTargChromName, StartLoci+1,MAPQ,szCigar,pszRNext,PNext,TLen);
NumEls = sscanf(pszSAMline,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%n",
         szDescriptor, &Flags, szChrom, &StartLoci,&MAPQ,szCigar,szRNext,&PNext,&TLen,&LineRemainIdx);
if(NumEls != 9)
	return(eBSFerrParse);

strncpy(pBAMalign->read_name,szDescriptor,sizeof(pBAMalign->read_name));	// szDescriptor as read name
pBAMalign->read_name[sizeof(pBAMalign->read_name)-1] = '\0';
pBAMalign->NumReadNameBytes = (INT32)strlen(pBAMalign->read_name) + 1;

pBAMalign->flag_nc = (UINT32)(Flags & 0x0ffff) << 16;						// flags

if(pBEDremapper != NULL)
	{
	if((ContigID = pBEDremapper->LocateFeatureIDbyName(szChrom)) < 1)
		return(eBSFerrFeature);
	pBEDremapper->GetFeature(ContigID,NULL,szChrom,&RelStartLoci);
	}
else
	RelStartLoci = 0;

strncpy(pBAMalign->szRefSeqName,szChrom,sizeof(pBAMalign->szRefSeqName));   // szChrom as szRefSeqName
pBAMalign->szRefSeqName[sizeof(pBAMalign->szRefSeqName)-1] = '\0';

if(m_szLastNotLocatedRefSeqName[0] != '\0')
	if(!stricmp(pBAMalign->szRefSeqName,m_szLastNotLocatedRefSeqName))
		return(eBSFerrFeature);

if((pBAMalign->refID = LocateRefSeqID(pBAMalign->szRefSeqName)) < 1)
	{
	strncpy(m_szLastNotLocatedRefSeqName,pBAMalign->szRefSeqName,sizeof(m_szLastNotLocatedRefSeqName) - 1);
	return(eBSFerrFeature);
	}

if(!(szRNext[0] == '*' || szRNext[0] == '=' ))
	{
	if(pBEDremapper != NULL)
		{
		if((ContigID = pBEDremapper->LocateFeatureIDbyName(szRNext)) < 1)
			return(eBSFerrFeature);
		pBEDremapper->GetFeature(ContigID,NULL,szRNext,&RelMateStartLoci);	
		if((pBAMalign->next_refID = LocateRefSeqID(szRNext)) < 1)
			return(eBSFerrFeature);
		}
	else
		RelMateStartLoci = 0;
	pBAMalign->next_pos = RelMateStartLoci + PNext - 1;
	pBAMalign->tlen = TLen;
	}
else
	{
	if(szRNext[0] == '=')
		{
		pBAMalign->next_refID = pBAMalign->refID;
		pBAMalign->next_pos = RelStartLoci + PNext - 1;
		pBAMalign->tlen = TLen;
		}
	else
		{
		pBAMalign->next_refID = -1;
		pBAMalign->next_pos = 0;
		pBAMalign->tlen = 0;
		}
	}
strncpy(pBAMalign->szMateRefSeqName,szRNext,sizeof(pBAMalign->szMateRefSeqName));   // szRNext as the mate SeqName
pBAMalign->szMateRefSeqName[sizeof(pBAMalign->szMateRefSeqName)-1] = '\0';

pBAMalign->pos = RelStartLoci + StartLoci - 1;								// StartLoci needs to be decremented
pBAMalign->end = pBAMalign->pos;

pBAMalign->bin_mq_nl = (MAPQ & 0x0ff) << 8;									// mapping quality
pBAMalign->bin_mq_nl |= pBAMalign->NumReadNameBytes & 0x0ff;

char *pCigar;
char CigarChr;
int CigarOP;
int CigarOPlen;
int CigarIdx;

CigarIdx = 0;
pCigar = szCigar;
CigarOPlen = 0;
if(*pCigar != '*')
	{
	while((CigarChr = *pCigar++) != '\0')
		{
		if(CigarIdx >= 50)
			return(eBSFerrParse);
		if(CigarChr >= '0' && CigarChr <= '9')
			{
			CigarOPlen *= 10;
			CigarOPlen += (int)(CigarChr - '0');
			continue;
			}
		switch(tolower(CigarChr)) {
			case 'm':
				CigarOP = 0;
				break;
			case 'i':
				CigarOP = 1;
				break;
			case 'd':
				CigarOP = 2;
				break;
			case 'n':
				CigarOP = 3;
				break;
			case 's':
				CigarOP = 4;
				break;
			case 'h':
				CigarOP = 5;
				break;
			case 'p':
				CigarOP = 6;
				break;
			case '=':
				CigarOP = 7;
				break;
			case 'x':				
				CigarOP = 8;
				break;
			default:				// unsupported!
				return(eBSFerrParse);
			}
		pBAMalign->cigar[CigarIdx++] = (UINT32)CigarOPlen << 4 | CigarOP;
		pBAMalign->end += CigarOPlen;
		if(CigarIdx > 50)
			return(eBSFerrParse);
		CigarOPlen = 0;
		CigarOP = 0;
		}
	pBAMalign->end -= 1;
	}
pBAMalign->flag_nc |= CigarIdx;
pBAMalign->NumCigarBytes = max(CigarIdx,1) * sizeof(UINT32);

pszSeq = &pszSAMline[LineRemainIdx];

char szSeq[cMaxReadLen+1];
char szQScore[cMaxReadLen+1];

sscanf(&pszSAMline[LineRemainIdx],"%s\t%s\t%n",szSeq,szQScore,&LineRemainIdx);
pBAMalign->l_seq = (INT32)strlen(szSeq);

if(szQScore[0] == '*')		// if there were no associated quality scores
	memset(pBAMalign->qual,0x0ff,pBAMalign->l_seq);
else
	strcpy((char *)pBAMalign->qual,szQScore);

UINT8 Byte;
int Ofs;

char Base;

Byte = 0;
Ofs = 0;
do
	{
	Base = szSeq[Ofs];
	switch(Base) {       // `=ACMGRSVTWYHKDBN' --> [0, 15];
		case 'a': case 'A':
			Byte |= 1;
			break;
		case 'c': case 'C':
			Byte |= 2;
			break;
		case 'g': case 'G':
			Byte |= 4;
			break;
		case 't': case 'T':
			Byte |= 8;
			break;
		default:
			Byte |= 15;
		}
	if(!(Ofs & 0x01))
		Byte <<= 4;

	if((Ofs & 0x01) || Ofs == pBAMalign->l_seq-1)
		{
		pBAMalign->seq[Ofs/2] = Byte;
		Byte = 0;
		}
	}
while(++Ofs < pBAMalign->l_seq);
pBAMalign->NumSeqBytes = (pBAMalign->l_seq + 1)/2;

return(eBSFSuccess);
}



// returns next SAM formated line
// if reading BAM input and special EOF block encountered then that block is simply sloughed as there could be blocks following
//


int											// number of chars returned in pszNxtLine
CSAMfile::GetNxtSAMline(char *pszNxtLine)	// copy next line read from input source to this line buffer; assumes caller has allocated at least cMaxBAMLineLen chars to hold worst case length
{
int Idx;
char Char;
int LenRead;
int LenRemaining;
int NxtLineLen;
UINT32 StartInBamIdx;	// m_pBAM[] at which current alignment entry starts 
UINT32 block_size;		// Length of the remainder of this alignment record (includes any auxiliary data)
INT32 refID;			// Reference sequence ID,  -1 <= refID < n_ref; -1 for a read without a mapping position
INT32 pos;				// 0-based leftmost coordinate (= POS - 1) 
UINT32 bin_mq_nl;		// bin<<16|MAPQ<<8|l_read_name ; bin is computed by the reg2bin(); l_read_name is the length of read name below (= length(QNAME) + 1).
UINT32 flag_nc;			// FLAG<<16|n_cigar_op; n_cigar_op is the number of operations in CIGAR
INT32 l_seq;			// Length of SEQ
INT32 next_refID;		// Ref-ID of the next segment (-1 <= mate_refID < n_ref)
INT32 next_pos;			// 0-based leftmost pos of the next segment (= PNEXT - 1)
INT32 tlen;				// Template length (= TLEN)
char szReadName[cMaxDescrIDLen+1];	// read_name  | char[l_read name] | NULL terminated (QNAME plus a tailing `\0')
UINT32 cigar;			// CIGAR: op_len<<4|op. `MIDNSHP=X'!`012345678'
int CigarOpLen;			// current CIGAR op len
char CigarOp;			// current CIGAR op
int CigarLen;			// cuurent length of constructed szCigar
char szCigar[cMaxBAMCigarOps * 5];		// CIGAR string
int BaseIdx;
char szSeq[cMaxBAMSeqLen];	// sequence
char szQual[cMaxBAMSeqLen];	// Phred base quality scores

char *pszTargSeq;
char *pszRNext;

if(pszNxtLine == NULL)
	return(eBSFerrParams);
*pszNxtLine = '\0';
LenRemaining = (int)(m_CurBAMLen - m_CurInBAMIdx);

if(m_bInEOF && LenRemaining == 0)	// if previously processed last line of input then return 0 to show no more lines to return
	return(0);

if(m_SAMFileType == eSFTSAM || m_SAMFileType == eSFTSAMgz)
	{
	LenRead = 0;
	while((LenRemaining = (int)(m_CurBAMLen - m_CurInBAMIdx)) > 0 || !m_bInEOF) 
		{
		// read next line(s) sloughing lines which are whitespace only or header lines starting with '@' as their first none-whitespace character
		if(!m_bInEOF && LenRemaining < (10 * cMaxBAMLineLen))		// try and buffer at least a few lines of source
			{
			if(LenRemaining && m_CurInBAMIdx > 0)
				{
				memmove(m_pBAM,&m_pBAM[m_CurInBAMIdx],LenRemaining);
				m_CurInBAMIdx = 0;
				m_CurBAMLen = LenRemaining;
				}
			if(m_SAMFileType == eSFTSAMgz)
				LenRemaining = gzread(m_gzInSAMfile,&m_pBAM[m_CurBAMLen],(int)(m_AllocBAMSize - m_CurBAMLen));
			else
				LenRemaining = read(m_hInSAMfile,&m_pBAM[m_CurBAMLen],(int)(m_AllocBAMSize - m_CurBAMLen));
			if(LenRemaining > 0)
				m_CurBAMLen += LenRemaining;
			else
				m_bInEOF = true;
			}
		Char = (char)m_pBAM[m_CurInBAMIdx++];
		m_TotInBAMProc += 1;
		if(LenRead == 0 && Char <= ' ')			// slough any leading control/whitespace chars
			continue;
		if(Char == '\n' || Char == '\r')
			{
			*pszNxtLine = '\0';
			return(LenRead);
			}
		*pszNxtLine++ = Char;
		LenRead += 1;
		}
	return(0);									// 0 to show no more lines to return
	}
else  // reading from BAM
	{
	LenRead = 0;
	while((LenRemaining = (int)(m_CurBAMLen - m_CurInBAMIdx)) > 0 || !m_bInEOF)
		{
		if(!m_bInEOF && LenRemaining < (10 * cMaxBAMLineLen))		// try and keep buffer full
			{
			if(LenRemaining && m_CurInBAMIdx > 0)
				{
				memmove(m_pBAM,&m_pBAM[m_CurInBAMIdx],LenRemaining);
				m_CurInBAMIdx = 0;
				m_CurBAMLen = LenRemaining;
				}
			if((LenRemaining = (int)bgzf_read(m_pInBGZF,&m_pBAM[m_CurBAMLen],(int)(m_AllocBAMSize - m_CurBAMLen))) < 0)		// will be -1 if errors ...
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetNxtSAMline: Not a BAM format file '%s'",m_szSAMfileName);
				Reset();
				return(eBSFerrFileAccess);
				}
			if(LenRemaining > 0)
				m_CurBAMLen += LenRemaining;
			else
				m_bInEOF = true;
			}

		// generally assured of at least cMaxBAMLineLen bytes from input loaded into m_pBAM, if less then must be near EOF
		if(m_TotInBAMProc == 4)				// load header length
			{
			m_InBAMHdrLen = *(int *)&m_pBAM[m_CurInBAMIdx];
			if(m_CurBAMLen < min(m_AllocBAMSize,(size_t)m_InBAMHdrLen + 8))
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetNxtSAMline: Unable to load header of expected length %d from BAM format file '%s'",m_InBAMHdrLen,m_szSAMfileName);
				Reset();
				return(eBSFerrFileAccess);
				}
			m_TotInBAMProc += 4;
			m_CurInBAMIdx += 4;
			}

		// if still within header then return each header line
		if(m_TotInBAMProc < ((size_t)m_InBAMHdrLen + 8))
			{
			do {
				if(m_CurInBAMIdx == m_CurBAMLen)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetNxtSAMline: Unable to complete loading headers of expected total length %d from BAM format file '%s'",m_InBAMHdrLen,m_szSAMfileName);
					Reset();
					return(eBSFerrFileAccess);
					}
				Char = (char)m_pBAM[m_CurInBAMIdx++];
				m_TotInBAMProc += 1;
				if(LenRead == 0 && Char <= ' ')			// slough any leading control/whitespace chars
					continue;
				if(Char == '\n' || Char == '\r')
					{
					*pszNxtLine = '\0';
					return(LenRead);
					}
				*pszNxtLine++ = Char;
				LenRead += 1;
				}
			while((m_TotInBAMProc < (8 + (size_t)m_InBAMHdrLen)) &&  LenRead < cMaxBAMLineLen-1);
			if(LenRead > 0)
				{
				*pszNxtLine = '\0';
				return(LenRead);
				}
			}

		// if still within target sequences then accumulate these until parsing alignment section
		if(m_TotInBAMProc == m_InBAMHdrLen + 8)
			{
			m_NumBAMSeqNames = *(int *)&m_pBAM[m_CurInBAMIdx];
			m_CurInBAMIdx += 4;
			m_TotInBAMProc += 4;
			m_NumRefSeqNames = 0;
			m_CurRefSeqsLen = 0;
			}
		if(m_NumRefSeqNames < m_NumBAMSeqNames)
			{
			if((m_CurRefSeqsLen + 1000) > m_AllocRefSeqsSize)	// ensure sufficent memory allocated to hold this new reference sequence
				{
				UINT8 *pTmp;
				pTmp = (UINT8 *)realloc(m_pRefSeqs,m_AllocRefSeqsSize + cAllocRefSeqSize);
				if(pTmp == NULL)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetNxtSAMline: unable to realloc memory for reference sequence names");
					Reset();
					return(eBSFerrMem);
					}
				m_pRefSeqs = (tsRefSeq *)pTmp;
				m_AllocRefSeqsSize += cAllocRefSeqSize;
				}
			m_NumRefSeqNames += 1;
			m_pCurRefSeq = (tsRefSeq *)((UINT8 *)m_pRefSeqs + m_CurRefSeqsLen);
			m_pCurRefSeq->SeqID = m_NumRefSeqNames;
			m_pCurRefSeq->SeqNameLen = -1 + *(int *)&m_pBAM[m_CurInBAMIdx];
			m_CurInBAMIdx += 4;
			strcpy(m_pCurRefSeq->szSeqName,(char *)&m_pBAM[m_CurInBAMIdx]);
			m_CurInBAMIdx += 1 + m_pCurRefSeq->SeqNameLen;
			m_pCurRefSeq->SeqLen = *(int *)&m_pBAM[m_CurInBAMIdx];
			m_CurInBAMIdx += 4;
			m_TotInBAMProc +=  m_pCurRefSeq->SeqNameLen + 9;
			m_CurRefSeqsLen += m_pCurRefSeq->SeqNameLen + sizeof(tsRefSeq);
			m_pCurRefSeq = NULL;
			continue;
			}

		// if within alignment section then can build alignment text
		if(m_pCurRefSeq == NULL)
			{
			m_pCurRefSeq = m_pRefSeqs;
			m_CurRefSeqNameID = 1;
			}

		block_size = *(UINT32 *)&m_pBAM[m_CurInBAMIdx];
		m_CurInBAMIdx += 4;
		m_TotInBAMProc += 4;
		StartInBamIdx = (UINT32)m_CurInBAMIdx;
		refID = *(int *)&m_pBAM[m_CurInBAMIdx];			// -1 if unaligned
		m_CurInBAMIdx += 4;

		if(refID == -1)
			pszTargSeq = (char *)"*";
		else
			{
			while(m_CurRefSeqNameID < m_NumBAMSeqNames && refID != m_pCurRefSeq->SeqID-1)
				{
				m_pCurRefSeq = (tsRefSeq *)((UINT8 *)m_pCurRefSeq + m_pCurRefSeq->SeqNameLen + sizeof(tsRefSeq));
				m_CurRefSeqNameID += 1;
				}
			if(refID != m_pCurRefSeq->SeqID-1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetNxtSAMline: unable to locate matching reference sequence name for %d",refID);
				Reset();
				return(eBSFerrEntry);
				}
			pszTargSeq = m_pCurRefSeq->szSeqName;
			}
		pos = *(int *)&m_pBAM[m_CurInBAMIdx];			// -1 if unaligned
		m_CurInBAMIdx += 4;
		bin_mq_nl = *(int *)&m_pBAM[m_CurInBAMIdx];		
		m_CurInBAMIdx += 4;
		flag_nc = *(int *)&m_pBAM[m_CurInBAMIdx];
		m_CurInBAMIdx += 4;
		l_seq = *(int *)&m_pBAM[m_CurInBAMIdx];
		m_CurInBAMIdx += 4;
		next_refID = *(int *)&m_pBAM[m_CurInBAMIdx];
		if(next_refID == -1)
			pszRNext = (char *)"*";
		else
			{
			if(next_refID == 0 || next_refID == refID)
				pszRNext = (char *)"=";
			else
				{
				tsRefSeq *pNxtRefSeq = m_pRefSeqs; 
				for(Idx = 0; Idx < (int)m_NumBAMSeqNames; Idx++)
					{
					if(next_refID == pNxtRefSeq->SeqID - 1)
						break;
					pNxtRefSeq = (tsRefSeq *)((UINT8 *)pNxtRefSeq + pNxtRefSeq->SeqNameLen + sizeof(tsRefSeq));
					}
				if(Idx == m_NumBAMSeqNames)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"GetNxtSAMline: unable to locate matching reference sequence name for %d",refID);
					Reset();
					return(eBSFerrEntry);
					}
				pszRNext = pNxtRefSeq->szSeqName;
				}
			}

		m_CurInBAMIdx += 4;
		next_pos = *(int *)&m_pBAM[m_CurInBAMIdx];
		m_CurInBAMIdx += 4;
		tlen = *(int *)&m_pBAM[m_CurInBAMIdx];
		m_CurInBAMIdx += 4;
		strncpy(szReadName,(char *)&m_pBAM[m_CurInBAMIdx],sizeof(szReadName));
		szReadName[sizeof(szReadName)-1] = '\0';
		m_CurInBAMIdx += bin_mq_nl & 0x0ff;

		if(refID == -1 || pos == -1)					// unaligned?
			{
			m_CurInBAMIdx += 4 * max(1,(flag_nc & 0x00ff));	// unaligned so slough any cigar ops
			strcpy(szCigar,"*");
			CigarLen = 1;
			}	
		else                                            // was aligned so build a CIGAR 
			{
			CigarLen = 0;
			for(UINT32 CigarIdx = 0; CigarIdx < (flag_nc & 0x00ff); CigarIdx++)
				{
				cigar = *(UINT32 *)&m_pBAM[m_CurInBAMIdx];
				m_CurInBAMIdx += 4;
				CigarOp = m_CigarOpsMap[cigar  & 0x00f];
				CigarOpLen = (cigar >> 4) & 0x0fffffff;
				CigarLen += sprintf((char *)&szCigar[CigarLen],"%d%c",CigarOpLen,CigarOp);
				}
			}

		// unpack the query sequence
		for(Idx = 0; Idx < l_seq; Idx++)
			{
			if(!(Idx & 0x01))
				BaseIdx = m_pBAM[m_CurInBAMIdx] >> 4;
			else
				BaseIdx = m_pBAM[m_CurInBAMIdx++];
			szSeq[Idx] = m_BasesMap[BaseIdx & 0x0f];
			}
		if(Idx & 0x01)
			m_CurInBAMIdx += 1;
		szSeq[l_seq] = '\0';

		// process Phred base quality scores
		if(m_pBAM[m_CurInBAMIdx] == 0x0ff)
			{
			szQual[0] = '*';
			szQual[1] = '\0';
			}
		else
			{
			memcpy(szQual,&m_pBAM[m_CurInBAMIdx],l_seq);
			szQual[l_seq] = '\0';
			}
		m_CurInBAMIdx += l_seq;
		if((m_CurInBAMIdx - StartInBamIdx) != block_size)	// if not exhusted current block then must be optional tags
			{
			// TODO: tag processing ....
			m_CurInBAMIdx =  StartInBamIdx + block_size;
			}
		m_TotInBAMProc += block_size;
		NxtLineLen = sprintf(pszNxtLine,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s",
						szReadName,(flag_nc >> 16) & 0x0ffff,pszTargSeq,pos+1,(bin_mq_nl >> 8) & 0x0ff,szCigar,pszRNext,next_pos+1,tlen,szSeq,szQual);
		return(NxtLineLen);
		}
	}
return(0);
}

int										// creat and initiate processing for SAM or BAM - with optional index - file generation
CSAMfile::Create(eSAMFileType SAMType,	// file type, expected to be either eSFTSAM or eSFTBAM_BAI or eSFTBAM_CSI 
				char *pszSAMFile,		// SAM(gz) or BAM file name
				int ComprLev,			// if BAM then BGZF compress at this requested level (0..9)
				char *pszVer)			// version text to use in generated SAM/BAM headers - if NULL then defaults to cszProgVer
{
if(SAMType < eSFTSAM || SAMType > eSFTBAM_CSI || pszSAMFile == NULL || pszSAMFile[0] == '\0')
	return(eBSFerrParams);

if((SAMType == eSFTSAM || SAMType >= eSFTBAM_BAI) && (ComprLev < 0 || ComprLev > 9))
	return(eBSFerrParams);

Reset();
if(pszVer != NULL && pszVer[0] != '\0')
	strncpy(m_szVer,pszVer,sizeof(m_szVer)-1);
else
	strcpy(m_szVer,cszProgVer);

m_SAMFileType = SAMType;
m_ComprLev = ComprLev;;
strcpy(m_szSAMfileName,pszSAMFile);
if(SAMType >= eSFTBAM_BAI)
	strcpy(m_szBAIfileName,pszSAMFile);

// try to create files
if(SAMType != eSFTSAMgz)
	{
#ifdef _WIN32
	m_hOutSAMfile = open(m_szSAMfileName,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutSAMfile = open(m_szSAMfileName,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutSAMfile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szSAMfileName,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hOutSAMfile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to create/truncate output file '%s'",m_szSAMfileName);
		Reset();
		return(eBSFerrCreateFile);
		}
	}
else
	{
	m_gzOutSAMfile = gzopen(m_szSAMfileName,"wb");
	if(m_gzOutSAMfile == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate gzip output file '%s'",m_szSAMfileName);
		return(eBSFerrCreateFile);
		}
	gzbuffer(m_gzOutSAMfile,cAllocSAMSize);		// large buffer to reduce number of writes required
	}


// files created/truncated
// make initial mem allocation
if(SAMType == eSFTSAM || SAMType == eSFTSAMgz)
	m_AllocBAMSize = cAllocSAMSize;
else
	m_AllocBAMSize = cAllocBAMSize;

if((m_pBAM = (UINT8 *)malloc(m_AllocBAMSize))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to alloc memory for bufferingr");
	Reset();
	return(eBSFerrMem);
	}

// if also BAI or CSI output then make initial mem allocations
if(SAMType >= eSFTBAM_BAI)
	{
	if((m_pBAI = (UINT8 *)malloc(cAllocBAISize))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to alloc memory for BAI");
		Reset();
		return(eBSFerrMem);
		}
	m_AllocBAISize = cAllocBAISize;

	if((m_pBAIChunks = (tsBAIChunk *)malloc(sizeof(tsBAIChunk) * cAllocBAIChunks))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to alloc memory for BAI or CSI chunks");
		Reset();
		return(eBSFerrMem);
		}
	m_AllocBAIChunks = cAllocBAIChunks;

	if((m_pChunkBins = (tsBAIbin *)malloc(sizeof(tsBAIbin) * cNumSAIBins))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to alloc memory for BAI or CSI chunks");
		Reset();
		return(eBSFerrMem);
		}
	memset(m_pChunkBins,0,sizeof(tsBAIbin) * cNumSAIBins);
	}

// alloc to hold reference sequence names
if((m_pRefSeqs = (tsRefSeq *)malloc(cAllocRefSeqSize))==NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to alloc memory for reference sequence names");
	Reset();
	return(eBSFerrMem);
	}
m_AllocRefSeqsSize = cAllocRefSeqSize;
m_NumRefSeqNames = 0;
m_CurRefSeqsLen = 0;
m_MaxRefSeqLen = 0;

m_CurBAILen = 0;
m_CurBAMLen = 0;

if(SAMType >= eSFTBAM)
	{
	// BAM will using BGZF compression ..
	char szLevel[10];
	sprintf(szLevel,"w%d",ComprLev);
	if((m_pBGZF = bgzf_dopen(m_hOutSAMfile, szLevel))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to initialise for BGZF compression (level %d) on file '%s'",m_szSAMfileName,ComprLev);
		Reset();
		return(eBSFerrMem);
		}
	m_hOutSAMfile = -1;
	m_pBAM[0] = (UINT8)'B';
	m_pBAM[1] = (UINT8)'A';
	m_pBAM[2] = (UINT8)'M';
	m_pBAM[3] = (UINT8)0x01;
	m_CurBAMLen = 8;				// the length of header text will be written into m_pBAM[4] when known
	}


m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"@HD\tVN:1.4\tSO:coordinate");

// initialisation completed
return(eBSFSuccess);
}

int					// returns current number of reference sequence names
CSAMfile::AddRefSeq(char *pszSpecies,	// sequence from this species
				  char *pszSeqName,		// sequence name
				  UINT32 SeqLen)		// sequence is of this length
{
tsRefSeq *pRefSeq;

if(m_SAMFileType >= eSFTBAM_BAI && (size_t)(UINT64)SeqLen > cMaxCSIRefSeqLen)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddREfSeq: Reference sequence '%s':'%s' of length %d is longer than max allowed %s index offset of %lld",
								pszSpecies,pszSeqName,SeqLen,
								m_SAMFileType == eSFTBAM_BAI ? "BAI" : "CSI",
								m_MaxIdxRefSeqLen-1);
	Reset();
	return(eBSFerrOfs);
	}

if((m_CurBAMLen + 1000) > m_AllocBAMSize)	// ensure sufficent memory allocated to hold this new reference sequence
	{
	UINT8 *pTmp;
	pTmp = (UINT8 *)realloc(m_pBAM,m_AllocBAMSize + cAllocBAMSize);
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to realloc memory for reference sequence names");
		Reset();
		return(eBSFerrMem);
		}
	m_pBAM = pTmp;
	m_AllocBAMSize += cAllocBAMSize;
	}

if((m_CurRefSeqsLen + 1000) > m_AllocRefSeqsSize)	// ensure sufficent memory allocated to hold this new reference sequence
	{
	UINT8 *pTmp;
	pTmp = (UINT8 *)realloc(m_pRefSeqs,m_AllocRefSeqsSize + cAllocRefSeqSize);
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to realloc memory for reference sequence names");
		Reset();
		return(eBSFerrMem);
		}
	m_pRefSeqs = (tsRefSeq *)pTmp;
	m_AllocRefSeqsSize += cAllocRefSeqSize;
	}

if(pszSpecies == NULL || pszSpecies[0] == '\0')
	m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\n@SQ\tAS:NA\tSN:%s\tLN:%u",pszSeqName,SeqLen);
else
	m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\n@SQ\tAS:%s\tSN:%s\tLN:%u",pszSpecies,pszSeqName,SeqLen);
m_NumBAMSeqNames += 1;

m_NumRefSeqNames += 1;
pRefSeq = (tsRefSeq *)((UINT8 *)m_pRefSeqs + m_CurRefSeqsLen);
pRefSeq->SeqID = m_NumRefSeqNames;
pRefSeq->SeqLen = SeqLen;
pRefSeq->SeqNameLen = (int)strlen(pszSeqName);
strcpy(pRefSeq->szSeqName,pszSeqName);
m_CurRefSeqsLen += sizeof(tsRefSeq) + pRefSeq->SeqNameLen;

if(SeqLen > m_MaxRefSeqLen)
	m_MaxRefSeqLen  = SeqLen;
return(m_NumBAMSeqNames);
}

int 
CSAMfile::StartAlignments(void)
{
UINT8 *pBAM;
tsRefSeq *pRefSeq;
size_t BGZFWritten;
UINT32 SeqNameID;

if(m_SAMFileType >= eSFTBAM_BAI)
	{
	if(m_SAMFileType == eSFTBAM_BAI && m_MaxRefSeqLen >= cMaxSAIRefSeqLen)
		{
		m_SAMFileType = eSFTBAM_CSI;
		gDiagnostics.DiagOut(eDLInfo,gszProcName,"StartAlignments: Generating CSI instead of SAI index file as alignments to sequence lengths (max %lld) more than SAI 512Mbp limit",m_MaxRefSeqLen);
		}
	if(m_SAMFileType == eSFTBAM_BAI)
		strcat(m_szBAIfileName,".bai");
	else
		strcat(m_szBAIfileName,".csi");
	}

if(m_SAMFileType >= eSFTBAM_BAI)
	{
	m_MaxIdxRefSeqLen = m_SAMFileType == eSFTBAM_BAI ? cMaxSAIRefSeqLen : cMaxCSIRefSeqLen;
	m_Alloc16KOfsVirtAddrsSize = m_MaxIdxRefSeqLen / 0x01000;
	
	if((m_p16KOfsVirtAddrs = (UINT64 *)malloc(m_Alloc16KOfsVirtAddrsSize))==NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"StartAlignments: unable to alloc %lld memory for BAI or CSI 16kb linear virtual offsets",m_Alloc16KOfsVirtAddrsSize);
		Reset();
		return(eBSFerrMem);
		}
	memset(m_p16KOfsVirtAddrs,0,m_Alloc16KOfsVirtAddrsSize);
	}

if(m_SAMFileType == eSFTBAM_BAI)
	{
#ifdef _WIN32
	m_hOutBAIfile = open(m_szBAIfileName,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
	if((m_hOutBAIfile = open(m_szBAIfileName,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOutBAIfile,0)!=0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szBAIfileName,strerror(errno));
			Reset();
			return(eBSFerrCreateFile);
			}
#endif
	if(m_hOutBAIfile < 0)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"StartAlignments: unable to create/truncate output file '%s'",m_szBAIfileName);
		Reset();
		return(eBSFerrCreateFile);
		}
	m_pBAI[0] = (UINT8)'B';
	m_pBAI[1] = (UINT8)'A';
	m_pBAI[2] = (UINT8)'I';
	m_pBAI[3] = (UINT8)0x01;
	m_CurBAILen = 4;
	}
else
	if(m_SAMFileType == eSFTBAM_CSI)
		{
		m_CSI_depth = CSIDepth(m_MaxRefSeqLen,m_CSI_min_shift);
		if(m_CSI_depth < 5)
			m_CSI_depth = 5;
#ifdef _WIN32
		m_hOutSAMfile = open(m_szBAIfileName,( O_WRONLY | _O_BINARY | _O_SEQUENTIAL | _O_CREAT | _O_TRUNC),(_S_IREAD | _S_IWRITE));
#else
		if((m_hOutSAMfile = open(m_szBAIfileName,O_WRONLY | O_CREAT,S_IREAD | S_IWRITE))!=-1)
			if(ftruncate(m_hOutSAMfile,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate %s - %s",m_szBAIfileName,strerror(errno));
				Reset();
				return(eBSFerrCreateFile);
				}
#endif
		if(m_hOutSAMfile < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"StartAlignments: unable to create/truncate output file '%s'",m_szBAIfileName);
			Reset();
			return(eBSFerrCreateFile);
			}
		char szLevel[10];
		sprintf(szLevel,"w%d",m_ComprLev);
		if((m_pgzOutCSIfile = bgzf_dopen(m_hOutSAMfile, szLevel))==NULL)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"StartAlignments: unable to initialise for BGZF compression (level %d) on file '%s'",m_szBAIfileName,m_ComprLev);
			Reset();
			return(eBSFerrMem);
			}
		m_hOutSAMfile = -1;
		m_pBAI[0] = (UINT8)'C';
		m_pBAI[1] = (UINT8)'S';
		m_pBAI[2] = (UINT8)'I';
		m_pBAI[3] = (UINT8)0x01;
		m_CurBAILen = 4;
		}

// if generating a BAM file then need to add the sequence names 
if(m_SAMFileType >= eSFTBAM)
	{
	m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\n@PG\tID:%s\tVN:%s",gszProcName,m_szVer); 
	*(UINT32 *)&m_pBAM[4] = (UINT32)(m_CurBAMLen - 8);
	*(UINT32 *)&m_pBAM[m_CurBAMLen] = m_NumRefSeqNames;
	m_CurBAMLen += 4;
	pRefSeq = m_pRefSeqs;
	for(SeqNameID = 0; SeqNameID < m_NumRefSeqNames; SeqNameID++)
		{
		if((m_CurBAMLen + 1000) > m_AllocBAMSize)	// ensure sufficent memory allocated to hold this new reference sequence
			{
			UINT8 *pTmp;
			pTmp = (UINT8 *)realloc(m_pBAM,m_AllocBAMSize + cAllocBAMSize);
			if(pTmp == NULL)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"StartAlignments: unable to realloc memory for reference sequence names");
				Reset();
				return(eBSFerrMem);
				}
			m_pBAM = pTmp;
			m_AllocBAMSize += cAllocBAMSize;
			}

		pBAM = &m_pBAM[m_CurBAMLen];
		*(UINT32 *)pBAM = pRefSeq->SeqNameLen + 1;
		pBAM += 4;
		strcpy((char *)pBAM,pRefSeq->szSeqName);
		pBAM += pRefSeq->SeqNameLen + 1;
		*(UINT32 *)pBAM = pRefSeq->SeqLen;
		m_CurBAMLen += 8 + pRefSeq->SeqNameLen + 1;
		pRefSeq = (tsRefSeq *)((INT8 *)pRefSeq + sizeof(tsRefSeq) + pRefSeq->SeqNameLen);
		}
	pBAM = &m_pBAM[m_CurBAMLen];
	if((BGZFWritten = (int)bgzf_write(m_pBGZF,m_pBAM,m_CurBAMLen))==-1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"StartAlignments: BGZF write failed");
		Reset();
		return(eBSFerrWrite);
		}
	m_CurBAMLen = 0;
	
	if(m_SAMFileType == eSFTBAM_BAI)
		{
		m_NumBAISeqNames = m_NumRefSeqNames;
		*(UINT32 *)&m_pBAI[m_CurBAILen] = m_NumBAISeqNames;
		m_CurBAILen += 4;
		}
	else
		if(m_SAMFileType == eSFTBAM_CSI)
			{
			*(UINT32 *)&m_pBAI[m_CurBAILen] = m_CSI_min_shift;
			m_CurBAILen += 4;
			*(UINT32 *)&m_pBAI[m_CurBAILen] = m_CSI_depth;
			m_CurBAILen += 4;
			*(UINT32 *)&m_pBAI[m_CurBAILen] = 0;
			m_CurBAILen += 4;
			m_NumBAISeqNames = m_NumRefSeqNames;
			*(UINT32 *)&m_pBAI[m_CurBAILen] = m_NumBAISeqNames;
			m_CurBAILen += 4;
			}
	}
else        // SAM output
	{
	m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\n@PG\tID:%s\tVN:%s\n",gszProcName,m_szVer);
	if(m_SAMFileType == eSFTSAM)
		CUtility::SafeWrite(m_hOutSAMfile,m_pBAM,m_CurBAMLen);
	else
		CUtility::SafeWrite_gz(m_gzOutSAMfile,m_pBAM,m_CurBAMLen);
	}
m_CurBAMLen = 0;
m_NumOf16Kbps = 0;
m_NumChunks = 0;
m_NumBinsWithChunks = 0;
m_pCurRefSeq = NULL;
return(m_NumBAMSeqNames);
}

// BAM index
// UINT8 magic[4];    // "BAI\1"
// UINT32 n_rf;       // number of reference sequences following
//    UINT32 n_bin;   // number of distinct bins for current reference sequence
//        UINT32 bin; // distinct bin
//        UINT32 chunks; // number of chunks following
//            UINT64 chumk_beg;		// virtual file offset at which chunk starts
//            UINT64 chumk_end;		// virtual file offset at which chunk ends
//    UINT32 n_intv;  // number of 16kb intervals for linear index
//        UINT64 ioffset;   // virtual file offset of first alignment in interval


// CSI index
// UINT8 magic[4];    // "CSI\1"
// UINT32 min_shift;  // # bits for minimum interval - defaults to 14 which is same as used in the BAI indexes
// UINT32 depth;      // R-tree depth - defaults to 5 which is same as used in the BAI indexes
// UINT32 l_aux;      // auxilary data - not currently used so defaults to 0
// UINT8  aux[l_aux]; // not currently required as no auxilary data used
// UINT32 n_rf;       // number of reference sequences following, same as in BAI
//    UINT32 n_bin;   // number of distinct bins for current reference sequence
//        UINT32 bin; // distinct bin
//        UINT64 loffset;   // virtual file offset of first overlapping record for bin
//        UINT32 chunks; // number of chunks following
//            UINT64 chumk_beg;		// virtual file offset at which chunk starts
//            UINT64 chumk_end;		// virtual file offset at which chunk ends

int									// write index to disk, returns number of bytes written, can be 0 if none attempted to be written, < 0 if errors
CSAMfile::WriteIdxToDisk(void)	
{
int BGZFWritten = (int)m_CurBAILen;
if(m_CurBAILen == 0 || m_SAMFileType < eSFTBAM_BAI)
	return(0);
if((m_hOutBAIfile == -1 && m_SAMFileType == eSFTBAM_BAI) || (m_pgzOutCSIfile == NULL && m_SAMFileType == eSFTBAM_CSI) || m_pBAI == NULL)
	return(eBSFerrFileClosed);

if(m_SAMFileType == eSFTBAM_BAI)
	{
	if(!CUtility::SafeWrite(m_hOutBAIfile,m_pBAI,m_CurBAILen))
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"UpdateSAIIndex: write to '%s' failed",m_szBAIfileName);
		Reset();
		return(eBSFerrWrite);
		}
	}
else
	{
	if((BGZFWritten = (int)bgzf_write(m_pgzOutCSIfile,m_pBAI,m_CurBAILen))==-1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"UpdateSAIIndex: write to '%s' failed",m_szBAIfileName);
		Reset();
		return(eBSFerrWrite);
		}
	}
m_CurBAILen = 0;
return(BGZFWritten);
}


int
CSAMfile::UpdateSAIIndex(bool bFinal)	// true if this is the final index update
{
int Rslt;
int BinIdx;
int ChunkIdx;
UINT32 *pSAI;
tsBAIbin *pBAIbin;
tsBAIChunk *pBAIChunks;

if((m_CurBAILen + 1000) > m_AllocBAISize)
	{
	if((Rslt = WriteIdxToDisk()) < eBSFSuccess)
		return(Rslt);
	m_CurBAILen = 0;
	}

pSAI = (UINT32 *)&m_pBAI[m_CurBAILen];

*pSAI++ = m_NumBinsWithChunks; // n_bin
m_CurBAILen += 4;

if(m_NumBinsWithChunks > 0)
	{
	pBAIbin = m_pChunkBins;
	for(BinIdx = 0; BinIdx < cNumSAIBins; BinIdx++,pBAIbin++)
		{
		if(pBAIbin->NumChunks)
			{
			*pSAI++ = BinIdx;	// distinct bin
			m_CurBAILen += 4;
			pBAIChunks = &m_pBAIChunks[pBAIbin->FirstChunk];
			if(m_SAMFileType == eSFTBAM_CSI)
				{
				*(UINT64 *)pSAI = pBAIbin->StartVA;
				pSAI += 2;
				m_CurBAILen += 8;
				}
			*pSAI++ = pBAIbin->NumChunks; // number of chunks following
			m_CurBAILen += 4;
			for(ChunkIdx =0;ChunkIdx < (int)pBAIbin->NumChunks;ChunkIdx++)
				{
				*(UINT64 *)pSAI = pBAIChunks->StartVA;
				pSAI += 2;
				*(UINT64 *)pSAI = pBAIChunks->EndVA;
				pSAI += 2;
				m_CurBAILen += 16;
				pBAIChunks = &m_pBAIChunks[pBAIChunks->NextChunk];

				if((m_CurBAILen + 1000) > m_AllocBAISize)
					{
					if((Rslt = WriteIdxToDisk()) < eBSFSuccess)
						return(Rslt);
					m_CurBAILen = 0;
					pSAI = (UINT32 *)m_pBAI; 
					}
				}
			}
		}

	if((m_CurBAILen + 1000 + (m_NumOf16Kbps * sizeof(UINT64))) > m_AllocBAISize)
		{
		if((Rslt = WriteIdxToDisk()) < eBSFSuccess)
			return(Rslt);
		m_CurBAILen = 0;
		pSAI = (UINT32 *)m_pBAI; 
		}
	if(m_SAMFileType == eSFTBAM_BAI)
		{
		*pSAI++ = m_NumOf16Kbps;
		m_CurBAILen += 4;
		memcpy(pSAI,m_p16KOfsVirtAddrs,m_NumOf16Kbps * sizeof(UINT64));
		pSAI += m_NumOf16Kbps * 2;
		m_CurBAILen += m_NumOf16Kbps * sizeof(UINT64);
		}
	}
else     // no alignments to this sequence
	{
	*pSAI = 0; // n_intv
	m_CurBAILen += 4;	
	}

if((m_CurBAILen + 1000) > m_AllocBAISize)
	{
	if((Rslt = WriteIdxToDisk()) < eBSFSuccess)
		return(Rslt);
	m_CurBAILen = 0;
	}
return(eBSFSuccess);
}

int
CSAMfile::AddChunk(UINT64 StartVA,		// start alignment BAM record is at this virtual address
				UINT32 Start,			// chunk starts at this loci
				UINT64 EndVA,			// chunk alignment BAM record ends at this virtual address
				UINT32 End)				// chunk ends at this loci
{
int Bin;
int KOfs;
tsBAIChunk *pBAIChunk;
tsBAIbin *pBin;

if(End >= m_MaxIdxRefSeqLen)		// exceeding the limit should have been picked up earlier but better safe than sorry
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddAlignment: Unable to accept alignment ending at %u in %s index extending past %lldbp'",
									End,m_SAMFileType == eSFTBAM_BAI ? "BAI" : "CSI", m_MaxIdxRefSeqLen);
	Reset();
	return(eBSFerrOfs);
	}

// if could potentially need more chunks than currently allocated then allocate
if((m_NumChunks + 5) >= m_AllocBAIChunks)	
	{
	tsBAIChunk *pTmp;
	pTmp = (tsBAIChunk *)realloc(m_pBAIChunks,sizeof(tsBAIChunk) * (m_AllocBAIChunks + cAllocBAIChunks));
	if(pTmp == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Create: unable to realloc memory for BAI chunks");
		Reset();
		return(eBSFerrMem);
		}
	m_pBAIChunks = pTmp;
	m_AllocBAIChunks += cAllocBAIChunks;
	}

KOfs = Start/0x04000;
if(m_p16KOfsVirtAddrs[KOfs] == 0)
	{
	m_NumOf16Kbps = KOfs + 1;
	m_p16KOfsVirtAddrs[KOfs] = StartVA;
	}

if(m_SAMFileType == eSFTBAM_BAI)
	Bin = BAIreg2bin(Start,End);	// which bin contains this chunk?
else
	Bin = CSIreg2bin(Start,End,m_CSI_min_shift,m_CSI_depth);	// which bin contains this chunk?	
pBin = &m_pChunkBins[Bin];
if(pBin->NumChunks == 0)		// first chunk allocated for this bin?
	{
	m_NumBinsWithChunks += 1;
	pBAIChunk = &m_pBAIChunks[m_NumChunks++];
	pBin->NumChunks = 1;
	pBin->StartVA = StartVA;
	pBin->FirstChunk = m_NumChunks-1;
	pBin->LastChunk = m_NumChunks-1;
	pBAIChunk->Bin = Bin;
	pBAIChunk->NextChunk = 0;	// first and currently only chunk in this bin
	pBAIChunk->Start = Start;
	pBAIChunk->StartVA = StartVA;
	pBAIChunk->End = End;
	pBAIChunk->EndVA = EndVA;
	return(m_NumChunks);
	}

// must have already been at least one chunk allocated to bin
if(pBin->StartVA > StartVA)
	pBin->StartVA = StartVA;

pBAIChunk = &m_pBAIChunks[pBin->LastChunk];
if(Start > (pBAIChunk->End + 1))
	{
	// starting a new chunk in this bin...
	pBAIChunk->NextChunk = m_NumChunks;
	pBAIChunk = &m_pBAIChunks[m_NumChunks++];
	pBin->NumChunks += 1;
	pBin->LastChunk = m_NumChunks-1;
	pBAIChunk->Bin = Bin;
	pBAIChunk->NextChunk = 0;
	pBAIChunk->Start = Start;
	pBAIChunk->StartVA = StartVA;
	pBAIChunk->End = End;
	pBAIChunk->EndVA = EndVA;
	}
else
	{
	// continue with same chunk, although chunk may need extending
	if(pBAIChunk->Start > Start)
		{
		pBAIChunk->Start = Start;
		pBAIChunk->StartVA = StartVA;
		}
	if(pBAIChunk->End < End)
		pBAIChunk->End = End;
	pBAIChunk->EndVA = EndVA;
	}
return(m_NumChunks);
}


// following BAM bin functions are copied from the specification at http://samtools.sourceforge.net/SAMv1.pdf
/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
int CSAMfile::BAIreg2bin(int beg, int end)
{
--end;
if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
#define MAX_BIN (((1<<18)-1)/7)
int CSAMfile::BAIreg2bins(int beg, int end, UINT16 *plist)
{
int i = 0, k;
--end;
plist[i++] = 0;
for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) plist[i++] = k;
for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) plist[i++] = k;
for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) plist[i++] = k;
for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) plist[i++] = k;
for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) plist[i++] = k;
return i;
}


// following BAM bin functions are copied from the specification at http://samtools.sourceforge.net/CSIv1.pdf
/* calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open) */
int CSAMfile::CSIreg2bin(INT64 beg,     // begins at inclusive of start 
			INT64 end,					// ends at 
			 int min_shift,				// # bits for minimum interval - defaults to 14 which is same as used in the BAI indexes
			 int depth)					// R-tree depth - defaults to 5 which is same as used in the BAI indexes
{
int l, s = min_shift, t = ((1<<depth*3) - 1) / 7;
for (--end, l = depth; l > 0; --l, s += 3, t -= 1<<l*3)
	if (beg>>s == end>>s) 
		return t + (int)((beg>>s));
return 0;
}

/* calculate the list of bins that may overlap with region [beg,end) (zero-based) */
int CSAMfile::CSIreg2bins(INT64 beg,		// begins at inclusive of start 
			INT64 end,						// ends at 
			int min_shift,					// # bits for minimum interval - defaults to 14 which is same as used in the BAI indexes
			int depth,						// depth (R-tree levels) from the maximum targeted sequence length
			int *bins)						// where to retun bins
{
int l, t, n, s = min_shift + depth*3;
for (--end, l = n = t = 0; l <= depth; s -= 3, t += 1<<l*3, ++l) 
	{
	int b = t + (int)(beg>>s), e = t + (int)(end>>s), i;
	for (i = b; i <= e; ++i) 
		bins[n++] = i;
	}
return n;
}

/* calulate depth (R-tree levels) from the maximum targeted sequence length */
int CSAMfile::CSIDepth(INT64 MaxSeqLen,			// max sequence length for any alignment ending at 
			 int min_shift)						// # bits for minimum interval - defaults to 14 which is same as used in the BAI indexes
{
INT64 s;
int n_lvls;
for (n_lvls = 0, s = (INT64)(1 << min_shift); MaxSeqLen > s; ++n_lvls, s <<= 3);
return(n_lvls);
}


int
CSAMfile::AddAlignment(tsBAMalign *pBAMalign,   // alignment to report
		  bool bLastAligned)					// true if this is the last read which was aligned, may be more reads but these are non-aligned reads
{
int Rslt;
int BGZFWritten;
UINT8 AlignBlock[10000];		// hold a single alignment!
UINT8 *pAlignBlock;
UINT64 CurStartVirtAddress;
UINT64 CurEndVirtAddress;

if(m_SAMFileType == eSFTSAM || m_SAMFileType == eSFTSAMgz)
	{
	char CigarOp;
	int OpLen;
	if(m_CurBAMLen + (cMaxFastQSeqLen * 10) > m_AllocBAMSize)
		{
		if(m_SAMFileType == eSFTSAM)
			CUtility::SafeWrite(m_hOutSAMfile,m_pBAM,m_CurBAMLen);
		else
			CUtility::SafeWrite_gz(m_gzOutSAMfile,m_pBAM,m_CurBAMLen);
		m_CurBAMLen = 0;
		}

	m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%s\t%d\t%s\t%d\t%d\t",
				pBAMalign->read_name,(pBAMalign->flag_nc >> 16) & 0x00ffff,pBAMalign->szRefSeqName,pBAMalign->pos+1,(pBAMalign->bin_mq_nl >> 8) & 0x0ff);

	if((pBAMalign->flag_nc & 0x00ff) == 0)
		*(char *)&m_pBAM[m_CurBAMLen] = '*';
	else
		{
		for(UINT32 CigarIdx = 0; CigarIdx < (pBAMalign->flag_nc & 0x00ff); CigarIdx++)
			{
			CigarOp = m_CigarOpsMap[pBAMalign->cigar[CigarIdx]  & 0x00f];
			OpLen = (pBAMalign->cigar[CigarIdx] >> 4) & 0x0fffffff;
			m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%d%c",OpLen,CigarOp);
			}
		}

	m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%c\t%d\t%d\t",
			pBAMalign->next_refID == -1 ? '*' : '=',pBAMalign->next_refID == -1 ? 0 : pBAMalign->next_pos + 1,pBAMalign->tlen);


	for(int SeqIdx = 0; SeqIdx < pBAMalign->l_seq; SeqIdx+=1)
		{
		if(!(SeqIdx & 0x01))
			m_pBAM[m_CurBAMLen++] = m_BasesMap[(pBAMalign->seq[SeqIdx/2] >> 4) & 0x0f];
		else
			m_pBAM[m_CurBAMLen++] = m_BasesMap[pBAMalign->seq[SeqIdx/2] & 0x0f];
		}
	m_pBAM[m_CurBAMLen++] = '\t';
	if(pBAMalign->qual[0] == 0xff)
		m_pBAM[m_CurBAMLen++] = '*';
	else
		{
		memcpy(&m_pBAM[m_CurBAMLen],pBAMalign->qual,pBAMalign->l_seq);
		m_CurBAMLen += pBAMalign->l_seq;
		}

	// optional alignment tags
	if(pBAMalign->NumAux > 0)
		{
		// current user tags supported are those for specifying reasons as to why a read was not accepted as being aligned
		// these tags are of the form:
		// YU:Z:<reason> where reason is the NAR enumeration as text
		int Idx;
		int Ofs;
		UINT8 *pVal;
		tsBAMauxData *pAuxData;
		pAuxData = pBAMalign->auxData;
		for(Idx = 0; Idx < pBAMalign->NumAux; Idx++,pAuxData++)
			{
			m_pBAM[m_CurBAMLen++] = '\t';
			switch(pAuxData->val_type) {    //  one of Value type: AcCsSiIfZHB
				case 'A':					//  type is a single printable char
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:A:%c",pAuxData->tag,pAuxData->value[0]);
					break;
				case 'B':					// type is an array of numerics which could be int8, uint8, int16, uint16, int32, uint32, or float
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:B:%c",pAuxData->tag,pAuxData->array_type);
					switch(pAuxData->array_type) {
						case 'c': 
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%1d",*(INT8 *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;
						case 'C':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%1u",*(UINT8 *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;
						case 's':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=2)
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%1d",*(INT16 *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;
						case 'S':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=2)
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%1u",*(UINT16 *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;
						case 'i':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=4)
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%1d",*(INT32 *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;
						case 'I':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=4)
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%1u",*(UINT32 *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;
						case 'f':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=sizeof(float))
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%f",*(float *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;
						default:
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
								{
								m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"%1d",*(INT8 *)pVal);
								if(Ofs < (pAuxData->NumVals - 1))
									m_pBAM[m_CurBAMLen++] = ',';
								}
							break;						
						}

				case 'H':					// type is a hex array
					for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
						{
						m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"0x%2x",*(UINT8 *)pVal);
						if(Ofs < (pAuxData->NumVals - 1))
							m_pBAM[m_CurBAMLen++] = ',';
						}
					break;
				case 'Z':					// printable string which may contain spaces
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:Z:%s",pAuxData->tag,(char *)pAuxData->value);
					break;
				case 'f':					// single float
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:f:%f",pAuxData->tag,*(float *)pAuxData->value);
					break;

				
				case 'i':					// INT32
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:i:%d",pAuxData->tag,*(INT32 *)pAuxData->value);
					break;
				case 'I':					// UINT32
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:i:%u",pAuxData->tag,*(UINT32 *)pAuxData->value);
					break;
				case 's':					// INT16
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:i:%d",pAuxData->tag,*(INT16 *)pAuxData->value);
					break;
				case 'S':					// UINT16
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:i:%u",pAuxData->tag,*(UINT16 *)pAuxData->value);
					break;
				case 'c':					// INT8
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:i:%d",pAuxData->tag,*(INT8 *)pAuxData->value);
					break;
				case 'C':					// UINT8
					m_CurBAMLen += sprintf((char *)&m_pBAM[m_CurBAMLen],"\t%.2s:i:%u",pAuxData->tag,*(UINT8 *)pAuxData->value);
					break;
				}
			}
		}
	m_pBAM[m_CurBAMLen++] = '\n';
	}
else   // BAM processing
	{
	// if read has been aligned then check that the alignment is within the limits for SAI or CSI
	if(pBAMalign->pos >= 0 && m_SAMFileType >= eSFTBAM_BAI && pBAMalign->end >= m_MaxIdxRefSeqLen)		
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddAlignment: Unable to accept alignment ending at %d in %s index extending past %lldbp",
			pBAMalign->end,m_SAMFileType == eSFTBAM_BAI ? "BAI" : "CSI",m_MaxIdxRefSeqLen-1);
		Reset();
		return(eBSFerrOfs);
		}

	// iterate over the reference sequence dictionary, from the current, looking for a match on the reference sequence
	if(pBAMalign->pos >= 0)
		{
		if(m_pCurRefSeq == NULL)	// NULL if first time
			{
			m_pCurRefSeq = m_pRefSeqs;
			if(m_SAMFileType >= eSFTBAM_BAI)
				{
				if(m_NumChunks)
					{
					memset(m_pChunkBins,0,sizeof(tsBAIbin) * cNumSAIBins);
					m_NumChunks = 0;
					}

				if(m_NumOf16Kbps)
					{
					memset(m_p16KOfsVirtAddrs,0,m_NumOf16Kbps * sizeof(UINT64));
					m_NumOf16Kbps = 0;
					}
				m_NumBinsWithChunks = 0;
				}
			}
		while(stricmp((char *)m_pCurRefSeq->szSeqName,pBAMalign->szRefSeqName))
			{
			if(m_pCurRefSeq->SeqID == m_NumRefSeqNames)		// unable to locate a match?
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddAlignment: Unable to locate matching reference sequence '%s'",pBAMalign->szRefSeqName);
				Reset();
				return(eBSFerrInternal);
				}
			// starting onto a new reference
			if(m_SAMFileType >= eSFTBAM_BAI)
				{
					// update index for previous
				UpdateSAIIndex();
				if(m_NumChunks)
					{
					memset(m_pChunkBins,0,sizeof(tsBAIbin) * cNumSAIBins);
					m_NumChunks = 0;
					}
				if(m_NumOf16Kbps)
					{
					memset(m_p16KOfsVirtAddrs,0,m_Alloc16KOfsVirtAddrsSize);
					m_NumOf16Kbps = 0;
					}
				m_NumBinsWithChunks = 0;
				}
			m_pCurRefSeq = (tsRefSeq *)((UINT8 *)m_pCurRefSeq + m_pCurRefSeq->SeqNameLen + sizeof(tsRefSeq));
			}
		pBAMalign->refID = m_pCurRefSeq->SeqID - 1;
		}
	else   // non-aligned read
		pBAMalign->refID = -1;

	pAlignBlock = &AlignBlock[4];
	*(UINT32 *)pAlignBlock = pBAMalign->refID;
	pAlignBlock += 4;
	*(UINT32 *)pAlignBlock = pBAMalign->pos;
	pAlignBlock += 4;
	*(UINT32 *)pAlignBlock =pBAMalign->bin_mq_nl;
	pAlignBlock += 4;
	*(UINT32 *)pAlignBlock =pBAMalign->flag_nc;
	pAlignBlock += 4;
	*(UINT32 *)pAlignBlock =pBAMalign->l_seq;
	pAlignBlock += 4;
	*(UINT32 *)pAlignBlock = pBAMalign->next_refID == -1 ? -1 : pBAMalign->refID;
	pAlignBlock += 4;
	*(UINT32 *)pAlignBlock =pBAMalign->next_pos;
	pAlignBlock += 4;
	*(UINT32 *)pAlignBlock =pBAMalign->tlen;
	pAlignBlock += 4;

	memcpy(pAlignBlock,pBAMalign->read_name,pBAMalign->NumReadNameBytes);
	pAlignBlock += pBAMalign->NumReadNameBytes; 

	memcpy(pAlignBlock,pBAMalign->cigar,pBAMalign->NumCigarBytes);
	pAlignBlock += pBAMalign->NumCigarBytes;
	memcpy(pAlignBlock,pBAMalign->seq,pBAMalign->NumSeqBytes);
	pAlignBlock += pBAMalign->NumSeqBytes;
	memcpy(pAlignBlock,pBAMalign->qual,pBAMalign->l_seq);
	pAlignBlock += pBAMalign->l_seq;

	// optional alignment tag(s) processing
	UINT32 AuxDataLen;
	AuxDataLen = 0;
	if(pBAMalign->NumAux > 0)
		{
		// user tags supported are:
		//	YU:Z:ML - the read exceeded the multiloci limit (ML) specified by the -R option
		//	YU:Z:NVA - the read had no valid alignments (NVA)
		int Idx;
		int Ofs;
		UINT8 *pVal;
		tsBAMauxData *pAuxData;
		pAuxData = pBAMalign->auxData;
		for(Idx = 0; Idx < pBAMalign->NumAux; Idx++,pAuxData++)
			{
			*pAlignBlock++ = pAuxData->tag[0];
			*pAlignBlock++ = pAuxData->tag[1];
			if(pAuxData->val_type == 'H')      // type is actually hex array - the SAM/BAM specification doesn't cover how hex arrays are to be formated when generating BAM so treating as if BYTE UINT8 array!
				*pAlignBlock++ = 'B';
			else
				*pAlignBlock++ = pAuxData->val_type;
			AuxDataLen += 3;

			switch(pAuxData->val_type) {    //  one of Value type: AcCsSiIfZHB
				case 'A':					//  type is a single printable char
					*pAlignBlock++ = pAuxData->value[0];
					AuxDataLen += 1;
					break;

				case 'B':					// type is an array of numerics which could be int8, uint8, int16, uint16, int32, uint32, or float
					*pAlignBlock++ = pAuxData->array_type;
					*(int *)pAlignBlock = pAuxData->NumVals;
					pAlignBlock += 4;
					AuxDataLen += 5;

					switch(pAuxData->array_type) {
						case 'c': 
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
								{
								*pAlignBlock++ = *(INT8 *)pVal;
								AuxDataLen += 1;
								}
							break;
						case 'C':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
								{
								*pAlignBlock++ = *(UINT8 *)pVal;
								AuxDataLen += 1;
								}
							break;
						case 's':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=2)
								{
								*(INT16 *)pAlignBlock = *(INT16 *)pVal;
								pAlignBlock += 2;
								AuxDataLen += 2;
								}
							break;
						case 'S':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=2)
								{
								*(UINT16 *)pAlignBlock = *(UINT16 *)pVal;
								pAlignBlock += 2;
								AuxDataLen += 2;
								}
							break;
						case 'i':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=4)
								{
								*(INT32 *)pAlignBlock = *(INT32 *)pVal;
								pAlignBlock += 4;
								AuxDataLen += 4;
								}
							break;
						case 'I':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=4)
								{
								*(UINT32 *)pAlignBlock = *(UINT32 *)pVal;
								pAlignBlock += 4;
								AuxDataLen += 4;
								}
							break;
						case 'f':
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal+=sizeof(float))
								{
								*(float *)pAlignBlock = *(float *)pVal;
								pAlignBlock += sizeof(float);
								AuxDataLen += sizeof(float);
								}
							break;
						default:
							for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
								{
								*pAlignBlock++ = *(INT8 *)pVal;
								AuxDataLen += 1;
								}
							break;						
						}

				case 'H':					// type is hex array - the SAM/BAM specification doesn't cover how hex arrays are to be formated when generating BAM so treating as if UINT8 array!
					*pAlignBlock++ = 'C';
					*(int *)pAlignBlock = pAuxData->NumVals;
					pAlignBlock += 4;
					AuxDataLen += 5;
					for(Ofs = 0; Ofs < pAuxData->NumVals; Ofs++,pVal++)
						{
						*pAlignBlock++ = *(INT8 *)pVal;
						AuxDataLen += 1;
						}
					break;

				case 'Z':					// printable string which may contain spaces
					{
					int StrLen = (int)strlen((char *)pAuxData->value);
					strcpy((char *)pAlignBlock,(char *)pAuxData->value);
					pAlignBlock += 1 + StrLen;
					AuxDataLen += 1 + StrLen;
					}
					break;
				case 'f':					// single float
					*(float *)pAlignBlock = *(float *)pAuxData->value;
					pAlignBlock += sizeof(float);
					AuxDataLen += sizeof(float);
					break;

				case 'i':					// INT32
					*(INT32 *)pAlignBlock = *(INT32 *)pAuxData->value;
					pAlignBlock += sizeof(INT32);
					AuxDataLen += sizeof(INT32);
					break;
				case 'I':					// UINT32
					*(UINT32 *)pAlignBlock = *(UINT32 *)pAuxData->value;
					pAlignBlock += sizeof(UINT32);
					AuxDataLen += sizeof(UINT32);
					break;
				case 's':					// INT16
					*(INT16 *)pAlignBlock = *(INT16 *)pAuxData->value;
					pAlignBlock += sizeof(INT16);
					AuxDataLen += sizeof(INT16);
					break;
				case 'S':					// UINT16
					*(UINT16 *)pAlignBlock = *(UINT16 *)pAuxData->value;
					pAlignBlock += sizeof(UINT16);
					AuxDataLen += sizeof(UINT16);
					break;
				case 'c':					// INT8
					*(INT8 *)pAlignBlock = *(INT8 *)pAuxData->value;
					pAlignBlock += sizeof(INT8);
					AuxDataLen += sizeof(INT8);
					break;
				case 'C':					// UINT8
					*(UINT8 *)pAlignBlock = *(UINT8 *)pAuxData->value;
					pAlignBlock += sizeof(UINT8);
					AuxDataLen += sizeof(UINT8);
					break;
				}
			}
		}

	UINT32 AlignBlockSize = 8 * sizeof(UINT32) + AuxDataLen + pBAMalign->NumReadNameBytes + pBAMalign->NumCigarBytes + pBAMalign->NumSeqBytes + pBAMalign->l_seq;
	*(UINT32 *)&AlignBlock[0] = (UINT32)(pAlignBlock - &AlignBlock[0]) - 4;

	if(m_SAMFileType >= eSFTBAM_BAI && pBAMalign->refID >= 0)
		CurStartVirtAddress = bgzf_tell(m_pBGZF);

	BGZFWritten = (int)bgzf_write(m_pBGZF,AlignBlock,AlignBlockSize + sizeof(UINT32));
	if(m_SAMFileType >= eSFTBAM_BAI && bLastAligned && pBAMalign->refID >= 0)
		bgzf_flush(m_pBGZF);

	if(m_SAMFileType >= eSFTBAM_BAI && pBAMalign->refID >= 0)
		{
		CurEndVirtAddress = bgzf_tell(m_pBGZF);
		if((Rslt=AddChunk(CurStartVirtAddress,pBAMalign->pos,CurEndVirtAddress,pBAMalign->end))<eBSFSuccess)
			{
			Reset();
			return(Rslt);
			}
		}
	}
return(eBSFSuccess);
}

int
CSAMfile::Close(void)
{
int Rslt;
if(m_SAMFileType == eSFTSAM || m_SAMFileType == eSFTSAMgz)
	{
	if(m_SAMFileType == eSFTSAM)
		{
		if(m_hOutSAMfile != -1)
			{
			if(m_CurBAMLen)
				CUtility::SafeWrite(m_hOutSAMfile,m_pBAM,m_CurBAMLen);
#ifdef _WIN32
			_commit(m_hOutSAMfile);
#else
			fsync(m_hOutSAMfile);
#endif
			close(m_hOutSAMfile);
			m_hOutSAMfile = -1;
			}
		}
	else
		{
		if(m_gzOutSAMfile != NULL)
			{
			if(m_CurBAMLen)
				CUtility::SafeWrite_gz(m_gzOutSAMfile,m_pBAM,m_CurBAMLen);
			gzclose(m_gzOutSAMfile);
			m_gzOutSAMfile = NULL;
			}
		}
	m_CurBAMLen = 0;
	}
else
	{
	if(m_SAMFileType >= eSFTBAM_BAI)
		{
		if((Rslt=UpdateSAIIndex(true)) != eBSFSuccess)
			return(Rslt);

		WriteIdxToDisk();

		if(m_SAMFileType == eSFTBAM_BAI && m_hOutBAIfile != -1)
			{
#ifdef _WIN32
			_commit(m_hOutBAIfile);
#else
			fsync(m_hOutBAIfile);
#endif
			close(m_hOutBAIfile);
			m_hOutBAIfile = -1;
			}
		else
			{
			if(m_pBGZF != NULL)
				{
				bgzf_flush(m_pBGZF);
				bgzf_close(m_pBGZF);
				m_pBGZF = NULL;
				}

			bgzf_flush(m_pgzOutCSIfile);
			bgzf_close(m_pgzOutCSIfile);
			m_pgzOutCSIfile = NULL;
			}
		}
	}
return(eBSFSuccess);
}
