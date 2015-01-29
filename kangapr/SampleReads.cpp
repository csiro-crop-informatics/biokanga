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

#include "SampleReads.h"


CSampleReads::CSampleReads(void)
{
m_bInitialised = false;
Reset(false);
}


CSampleReads::~CSampleReads(void)
{
Reset(false);
}

void
CSampleReads::Reset(bool bSync)
{
if(m_bInitialised)
	{
	m_P1InFasta.Close();

	m_P2InFasta.Close();

	if(m_hOut5Reads != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOut5Reads);
#else
			fsync(m_hOut5Reads);
#endif
		close(m_hOut5Reads);
		m_hOut5Reads = -1;
		}

	if(m_hOut3Reads != -1)
		{
		if(bSync)
#ifdef _WIN32
			_commit(m_hOut3Reads);
#else
			fsync(m_hOut3Reads);
#endif
		close(m_hOut3Reads);
		m_hOut3Reads = -1;
		}
	}
else
	{
	m_hOut3Reads = -1;			// file handle for output 3' sampled reads
	m_hOut5Reads = -1;			// file handle for output 5' sampled reads
	}

m_szIn5ReadsFile[0] = 0;	// 5' reads are in this file
m_szIn3ReadsFile[0] = 0;	// 3' reads are in this file
m_szOut5File[0] = 0;		// write 5' sampled reads to this file
m_szOut3File[0] = 0;		// write 3' sampled reads to this file


m_bSOLiD = false;			// if true then input reads are expected to be in SOLiD colour space
m_bIsFastq = false;			// if basespace then input reads are fastq
m_bAppendOut = false;		// if true then append if output files exist, otherwise trunctate existing output files 
m_SampleOfs = 0;				// start at this read instance (1 = first read)
m_SampleIncr = 0;				// sample every Nth from SampleOfs
m_SampleMax = 0;				// sample at most this number of reads
m_bInitialised = true;
}

// OpenFiles
// Initialise and open files for processing
int
CSampleReads::OpenFiles(void)
{

// check if files specified for open/create/append...
if(!m_bInitialised || m_szIn5ReadsFile[0] == '\0' || (m_PMode == ePMPairedEnd && m_szIn3ReadsFile[0] == '\0'))
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"OpenFiles: unexpected state error");
	return(eBSFerrInternal);
	}





if(m_bAppendOut)
	{
#ifdef _WIN32
	m_hOut5Reads = open(m_szOut5File,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
	m_hOut5Reads = open(m_szOut5File,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
	if(m_hOut5Reads == -1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/open output 5' PE reads file '%s'",m_szOut5File);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	// seek to end of file ready for appending
	lseek(m_hOut5Reads,0,SEEK_END);

	if(m_PMode == ePMPairedEnd)
		{
		if(m_szOut3File[0] != '\0')
			{
#ifdef _WIN32
			m_hOut3Reads = open(m_szOut3File,( _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT),(_S_IREAD | _S_IWRITE));
#else
			m_hOut3Reads = open(m_szOut3File,O_RDWR | O_CREAT,S_IREAD | S_IWRITE);
#endif
			if(m_hOut3Reads == -1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/open output 3' PE reads file '%s'",m_szOut3File);
				Reset(false);
				return(eBSFerrCreateFile);
				}
			// seek to end of file ready for appending
			lseek(m_hOut3Reads,0,SEEK_END);
			}
		}
	}
else
	{
#ifdef _WIN32
	m_hOut5Reads = open(m_szOut5File,O_CREATETRUNC);
#else
	if((m_hOut5Reads = open(m_szOut5File,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
		if(ftruncate(m_hOut5Reads,0)!=0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate output %s 5' PE reads - %s",m_szOut5File,strerror(errno));
				Reset(false);
				return(eBSFerrCreateFile);
				}
#endif
	if(m_hOut5Reads == -1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output 5' PE microcontigs file '%s'",m_szOut5File);
		Reset(false);
		return(eBSFerrCreateFile);
		}
	// seek to end of file ready for appending
	lseek(m_hOut5Reads,0,SEEK_END);

	if(m_PMode == ePMPairedEnd)
		{
		if(m_szOut3File[0] != '\0')
			{
#ifdef _WIN32
			m_hOut3Reads = open(m_szOut3File,O_CREATETRUNC);
#else
			if((m_hOut3Reads = open(m_szOut3File,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
				if(ftruncate(m_hOut3Reads,0)!=0)
					{
					gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to truncate output %s 3' PE reads - %s",m_szOut3File,strerror(errno));
					Reset(false);
					return(eBSFerrCreateFile);
					}

#endif
			if(m_hOut3Reads == -1)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Process: unable to create/truncate output 3' PE file '%s'",m_szOut3File);
				Reset(false);
				return(eBSFerrCreateFile);
				}
			// seek to end of file ready for appending
			lseek(m_hOut3Reads,0,SEEK_END);
			}
		}
	}
return(eBSFSuccess);
}

int
CSampleReads::GenSamples(etPMode PMode,		// process mode, single or paired end
			  bool bAppendOut,				// if true then append if output files exist, otherwise trunctate existing output files 
			  int SampleOfs,				// start at this read instance (1 = first read)
			  int SampleIncr,				// sample every Nth from SampleOfs
			  int SampleMax,				// sample at most this number of reads	  
			  char *pszIn5ReadsFile,		// 5' reads are in this file
			  char *pszIn3ReadsFile,		// 3' reads are in this file
			  char *pszOut5ReadsFile,		// write sampled 5' reads to this file
			  char *pszOut3ReadsFile)		// write sampled 5' reads to this file
{
int Rslt;
int CurFastaSeqs;
int NumSamples;
int P1ReadLen;
int P2ReadLen;
UINT8 *pP1RawReadsBuff;
UINT8 *pP2RawReadsBuff;
UINT8 *pP1QBuff;
UINT8 *pP2QBuff;
UINT8 szP1DescrBuff[500];
UINT8 szP2DescrBuff[500];
int NumDescrReads;
int P1DescrLen;
int P2DescrLen;

Reset();
m_PMode = PMode;
m_bAppendOut = bAppendOut;
m_SampleOfs = SampleOfs;				// start at this read instance (1 = first read)
m_SampleIncr = SampleIncr;				// sample every Nth from SampleOfs
m_SampleMax = SampleMax;				// sample at most this number of reads	
strcpy(m_szIn5ReadsFile,pszIn5ReadsFile);
if(m_PMode == ePMPairedEnd)
	strcpy(m_szIn3ReadsFile,pszIn3ReadsFile);
else
	m_szIn3ReadsFile[0] = '\0';
strcpy(m_szOut5File,pszOut5ReadsFile);
if(m_PMode == ePMPairedEnd)
	strcpy(m_szOut3File,pszOut3ReadsFile);
else
	m_szOut3File[0] = '\0';
if((pP1RawReadsBuff = new UINT8 [cAllocRawSeqLen]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for raw reads buffering...");
	Reset(false);
	return(eBSFerrMem);
	}
if((pP1QBuff = new UINT8 [cAllocRawSeqLen]) == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for quality score buffering...");
	delete pP1RawReadsBuff;
	Reset(false);
	return(eBSFerrMem);
	}
if((Rslt=(teBSFrsltCodes)m_P1InFasta.Open(pszIn5ReadsFile,true))!=eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszIn5ReadsFile,m_P1InFasta.ErrText((teBSFrsltCodes)Rslt),m_P1InFasta.GetErrMsg());
	delete pP1RawReadsBuff;
	delete pP1QBuff;
	Reset(false);
	return(Rslt);
	}
if(m_PMode == ePMPairedEnd)
	{
	if((pP2RawReadsBuff = new UINT8 [cAllocRawSeqLen]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for raw reads buffering...");
		delete pP1RawReadsBuff;
		delete pP1QBuff;
		Reset(false);
		return(eBSFerrMem);
		}
	if((pP2QBuff = new UINT8 [cAllocRawSeqLen]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"Unable to allocate memory for quality score buffering...");
		delete pP1RawReadsBuff;
		delete pP1QBuff;
		delete pP2RawReadsBuff;
		Reset(false);
		return(eBSFerrMem);
		}
	if((Rslt=(teBSFrsltCodes)m_P2InFasta.Open(pszIn3ReadsFile,true))!=eBSFSuccess)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"LoadRawReads: Unable to open '%s' [%s] %s",pszIn3ReadsFile,m_P2InFasta.ErrText((teBSFrsltCodes)Rslt),m_P2InFasta.GetErrMsg());
		delete pP1RawReadsBuff;
		delete pP1QBuff;
		delete pP2RawReadsBuff;
		delete pP2QBuff;
		Reset(false);
		return(Rslt);
		}
	}
else
	{
	pP2RawReadsBuff = NULL;
	pP2QBuff = NULL;
	}
NumDescrReads = 0;
CurFastaSeqs = 0;
NumSamples = 0;
while((Rslt = (teBSFrsltCodes)(P1ReadLen = m_P1InFasta.ReadSequence(pP1RawReadsBuff,cAllocRawSeqLen,true,false))) > eBSFSuccess)
	{
	if(P1ReadLen == eBSFFastaDescr)		// just read a descriptor line which would be as expected for multifasta or fastq
		{
		NumDescrReads += 1;
		P1DescrLen = m_P1InFasta.ReadDescriptor((char *)szP1DescrBuff,sizeof(szP1DescrBuff)-1);
		szP1DescrBuff[sizeof(szP1DescrBuff)-1] = '\0'; 
		P1ReadLen = m_P1InFasta.ReadSequence(pP1RawReadsBuff,cAllocRawSeqLen);
		if(m_PMode == ePMPairedEnd)
			{
			Rslt = (teBSFrsltCodes)(P2ReadLen = m_P2InFasta.ReadSequence(pP2RawReadsBuff,cAllocRawSeqLen,true,false));
			if(Rslt <= 	eBSFSuccess)
				{
				while(m_P2InFasta.NumErrMsgs())
					gDiagnostics.DiagOut(eDLFatal,gszProcName,m_P2InFasta.GetErrMsg());
				delete pP1RawReadsBuff;
				delete pP2RawReadsBuff;
				m_P1InFasta.Close();
				m_P2InFasta.Close();
				}
			if(P2ReadLen != eBSFFastaDescr)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"raw sequence file '%s' processing error: %s ",pszIn3ReadsFile, 
														Rslt > 0 ? "over length read" : "not a multifasta short reads or fastq file");
				delete pP1RawReadsBuff;
				delete pP2RawReadsBuff;
				m_P1InFasta.Close();
				m_P2InFasta.Close();
				return(eBSFerrParse);
				}
			P2DescrLen = m_P2InFasta.ReadDescriptor((char *)szP2DescrBuff,sizeof(szP2DescrBuff)-1);
			szP2DescrBuff[sizeof(szP2DescrBuff)-1] = '\0'; 
			}
		if(P1ReadLen < 0)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed from '%s'",NumDescrReads-1,pszIn5ReadsFile);
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szP1DescrBuff);
			delete pP1RawReadsBuff;
			m_P1InFasta.Close();
			if(m_PMode == ePMPairedEnd)
				{
				delete pP2RawReadsBuff;
				m_P2InFasta.Close();
				}
			return(eBSFerrParse);
			}
		if(m_PMode == ePMPairedEnd)
			{
			P2ReadLen = m_P2InFasta.ReadSequence(pP2RawReadsBuff,cAllocRawSeqLen);
			if(P2ReadLen < 0)
				{
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Problem parsing sequence after %d reads parsed from '%s'",NumDescrReads-1,pszIn3ReadsFile);
				gDiagnostics.DiagOut(eDLFatal,gszProcName,"Last descriptor parsed: %s",szP2DescrBuff);
				delete pP1RawReadsBuff;
				delete pP2RawReadsBuff;
				m_P1InFasta.Close();
				m_P2InFasta.Close();
				return(eBSFerrParse);
				}
			}
		}
	if(CurFastaSeqs++ < SampleOfs)
		continue;
	if(SampleIncr > 1 && ((CurFastaSeqs - SampleOfs) % SampleIncr))
		continue;
	NumSamples += 1;
	if(SampleMax > 0 && (NumSamples >= SampleMax))
		break;
	}
#ifdef _WIN32
_commit(m_hOut5Reads);
#else
fsync(m_hOut5Reads);
#endif
close(m_hOut5Reads);
m_hOut5Reads = -1;
if(m_PMode == ePMPairedEnd)
	{
	#ifdef _WIN32
	_commit(m_hOut3Reads);
	#else
	fsync(m_hOut3Reads);
	#endif
	close(m_hOut3Reads);
	m_hOut3Reads = -1;
	}
delete pP1RawReadsBuff;
delete pP1QBuff;
m_P1InFasta.Close();
if(m_PMode == ePMPairedEnd)
	{
	delete pP2RawReadsBuff;
	delete pP2QBuff;
	m_P2InFasta.Close();
	}
return(Rslt);
}