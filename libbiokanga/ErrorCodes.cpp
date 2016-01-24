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
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif

CErrorCodes::CErrorCodes(void)
{
ClearErrs();
}

CErrorCodes::~CErrorCodes(void)
{
}

// clears all error messages
void 
CErrorCodes::ClearErrs(void)
{
m_NumMsgs = 0;
m_OldestMsg = 0;
for(int Idx = 0; Idx < cBSFMaxErrMsgs; Idx++)
	m_szErrMsgs[Idx][0] = '\0';
}

int 
CErrorCodes::NumErrMsgs(void)			// returns number of error messages
{
return(m_NumMsgs);
}

char *
CErrorCodes::GetErrMsg(void)			// pops and returns oldest error message
{
char *pszErr = (char *)"success";				// assume no errors
if(m_NumMsgs)
	{
	m_NumMsgs--;
	pszErr = m_szErrMsgs[m_OldestMsg++];
	if(m_OldestMsg >= cBSFMaxErrMsgs)
		m_OldestMsg = 0;
	}
return(pszErr);
}

int										// number of error messages 
CErrorCodes::AddErrMsg(const char *pszSource,	// used to identify message source
				  const char *pszFormat,...) // message text (sprintf format) plus any parameters
{
va_list Args;
int LineLen;
char szDiag[cBSFMaxErrMsgLen+1];
char szSource[cBSFMaxErrMsgSrc+1];
#ifdef _WIN32
struct _timeb timebuffer;
#else
struct timeb timebuffer;
#endif
char *timeline;
char *pszErrMsg;

pszErrMsg = m_szErrMsgs[(m_OldestMsg + m_NumMsgs)%cBSFMaxErrMsgs];
if(m_NumMsgs == cBSFMaxErrMsgs)
	{
	m_OldestMsg++;
	if(m_OldestMsg >= cBSFMaxErrMsgs)
		m_OldestMsg = 0;
	}
else
	m_NumMsgs++;

strncpy(szSource,pszSource,cBSFMaxErrMsgSrc);
szSource[cBSFMaxErrMsgSrc] = '\0';
#ifdef _WIN32
_ftime(&timebuffer);
#else
ftime(&timebuffer);
#endif
timeline = ctime(&timebuffer.time);
va_start(Args, pszFormat );
#ifdef _WIN32
_vsnprintf(szDiag,cBSFMaxErrMsgLen,pszFormat,Args);
#else
vsnprintf(szDiag,cBSFMaxErrMsgLen,pszFormat,Args);
#endif

szDiag[cBSFMaxErrMsgLen] = '\0';
LineLen = sprintf(pszErrMsg,"[%.15s.%03d %.4s] %s: %s\n",&timeline[4],(int)timebuffer.millitm, &timeline[20],szSource,szDiag);
return(m_NumMsgs);
}
				

const char *
CErrorCodes::ErrText(teBSFrsltCodes ErrCode)
{
switch(ErrCode) {
	case eBSFSuccess: return("No errors");
	case eBSFFastaDescr: return("fasta descriptor line avail"); 
	case eBSFerrParams:	return("parameter error");
	case eBSFerrFileClosed:	return("File must be opened before calling this method");
	case eBSFerrOfs:	return("offset specified past end of data");
	case eBSFerrCvrtType:	return("unable to convert to specified data type");
	case eBSFerrEntryCreate:return("Can't seal entry if none started"); 
	case eBSFerrMem:		return("unable to alloc memory"); 
	case eBSFerrNotBioseq:	return("file exists but not a bioseq file");
	case eBSFerrNotFasta:	return("file exists but not a fasta file");
	case eBSFerrFileType:	return("file is a bioseq file, but does not contain specified type");
	case eBSFerrFileDPT:	return("file is a bioseq file, but does not contain Data Points"); 
	case eBSFerrOpnFile:	return("unable to open file");
	case eBSFerrCreateFile:	return("unable to create file"); 
	case eBSFerrExecFile:	return("unable to execute external process");
	case eBSFerrClosed:		return("file is closed");
	case eBSFerrFileVer:	return("file version error");
	case eBSFerrFileAccess:	return("file access (seek/read/write) failed");
	case eBSFerrRead:		return("file not opened for read");
	case eBSFerrWrite:		return("file not opened for write");
	case eBSFerrTypeChange:	return("data type changed");
	case eBSFerrUserType:	return("can't curently handle user defined dataset types"); 
	case eBSFerrConfParam:	return("illegal structure parameter specified");
	case eBSFerrFastaDescr:	return("can't locate fasta descriptor");
	case eBSFerrFastaSymb:	return("illegal fasta sequence symbol");
	case eBSFerrFastqChr:   return("illegal fastq character");
	case eBSFerrFastqSeqID:	return("problem in fastq sequence identifier");
	case eBSFerrFastqSeq:	return("problem in fastq sequence");
	case eBSFerrFastqDescr:	return("problem in fastq descriptor");
	case eBSFerrFastqQScores: return("problem in fastq quality scores");
	case eBSFerrLocField:	return("Expected field not present");
	case eBSFerrStructParm:	return("Error processing structural parameters file");
	case eBSFerrStructStep:	return("Unable to determine conformation for step");
	case eBSFerrCentroidParam: return("Error processing centroids file");
	case eBSFerrNumSrcFiles:return("too many source files");
	case eBSFerrNoSrcFiles:	return("no source file");
	case eBSFerrNumAlignBlks:return("can't add any more alignment blocks"); 
	case eBSFerrAlignBlk:	return("can't locate alignment block or block not started");
	case eBSFerrObj:		return("unable to instantiate internal object"); 
	case eBSFerrRefDataset:	return("reference dataset error");
	case eBSFerrDataPtType:	return("unsupported data point type");
	case eBSFerrNoFeatures:	return("feature file is empty - no features");
	case eBSFerrNoEntries:	return("bioseq file is empty - no entries"); 
	case eBSFerrMaxDirEls:	return("hit max datasegs");
	case eBSFerrMaxChroms:	return("hit max chromosomes");
	case eBSFerrMaxDatasets:return("hit max datasets limit");
	case eBSFerrMaxFeatures:return("hit max features");
	case eBSFerrMaxEntries:	return("hit max entries");
	case eBSFerrFeature:	return("can't locate feature");
	case eBSFerrExon:		return("can't locate exon");
	case eBSFerrIntron:		return("can't locate intron");
	case eBSFerrEntry:		return("can't locate entry");
	case eBSFerrDataset:	return("can't locate dataset");
	case eBSFerrChrom:		return("can't locate chromosome");
	case eBSFerrSpecies:	return("can't locate species");
	case eBSFerrParse:		return("error whilst parsing file - unexpected format"); 
	case eBSFerrGene:		return("BEDfile file does not contain gene details"); 
	case eBSFerrProbMatrices: return("Error processing transitional probabilities matrices file");
	case eBSFErrBase:		return("Unrecognised base");
	case eBSFerrFileOpened:	return("File must be closed before calling this method");
	case eBSFerrFileName:	return("Unable to parse filename for ref+rel species names");
	case eBSFerrFieldCnt:	return("CSV number of fields parsed not number expected");
	case eBSFerrRowCnt:		return("CSV number of rows parsed not number expected");
	case eBSFerrFieldID:	return("CSV field identifier outside of range");
	case eBSFerrNumRange:	return("CSV field numeric conversion range error");
	case eBSFerrQuoteErrs:	return("CSV Additional chars followed end of quoted text");
	case eBSFerrQuoteIncomplete: return("CSV field unbalanced quotes");
	case eBSFerrRegion:		return("requested region not supported");

	case eBSFerrDupGOTerm:	return("Duplicate GO:Term processed");
	case eBSFerrGOID:		return("requested GO:Ident not located");
	case eBSFerrGOTagVal:	return("requested GO tag value not located");

	case eBSErrSession:		return("session requested not active");
	case cBSFSyncObjErr:	return("unable to synchronise access to object");

	case cBSFSocketErr:     return("socket level error");
	case cBSFNWSProtErr:    return("Services protocol error");
	case eBSFerrInternal:	return("!!Internal software processing error - bug!!"); 
	default:
		return("Unrecognised result code");
	}
}

