/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif
CCSV2BED::CCSV2BED(void)
{
m_pStepValues = NULL;
m_NumStepsAllocd = 0;
m_NumSteps = 0;
}

CCSV2BED::~CCSV2BED(void)
{
if(m_pStepValues != NULL)
	delete m_pStepValues;
}

int 
CCSV2BED::process(char *pszCSV,			// csv file containing structual values
				  char *pszBEDorWIG,	// name of BED or WIG file to create/append
				  char *pszTrackName,	// track name 
				  char *pszDescr,		// description
				  char *pszChrom,		// chromosome 
				  char *pValName,		// entry name
				  unsigned int NumClasses,	// number of classes into which scored values are to be segmented
				  unsigned int Windowsize,	// window size to use for data mean smoothing (1 if no smoothing)
				  unsigned int EntryFld,	// which field in the csv file contains the entry identifier (1..n)
				  unsigned int EntryID,		// entry identifer to match
				  unsigned int EntryOfsFld,	// which field in the csv file contains the entry offset base to left of step
				  unsigned int StartOfs,	// starting offset (base to left of last step)
				  unsigned int EndOfs,		// ending offset (base to left of last step)
  				  int Remap,				// remap offsets by Remap delta
				  unsigned int ValFld,		// which field contains the value of interest
  				  bool bAppend,			// true to append on to existing file otherwise (false) truncate existing file
				  bool bWIG)			// true if to output as WIG file, default is as BED
{
FILE *pCSVStream;
char szLineBuff[512];
char *pDst,*pSrc;
char *pEntryFld;
char *pEntryOfsFld;
char *pValFld;

char chr;
int LineNum;
int Entry;
int FldPsn;
unsigned int CurEntryOfs;
double StepValue;

if(EndOfs && StartOfs > EndOfs)
	return(eBSFerrParams);

if((pCSVStream = fopen(pszCSV,"r"))==NULL)
	{
	AddErrMsg("CCSV2BED::process","Unable to open %s - %s",pszCSV,strerror(errno));
	return(eBSFerrOpnFile);
	}

LineNum = 0;
while(fgets(szLineBuff,sizeof(szLineBuff),pCSVStream)!= NULL)
	{
	LineNum++;
	if(strlen(szLineBuff) < 1)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	// strip any whitespace and quotes
	pDst = pSrc = szLineBuff;
	while(chr = *pSrc++)
		if(!isspace(chr) && chr != '\'' && chr != '"')
			*pDst++ = chr;
	*pDst = '\0';
	if(szLineBuff[0] == '\0')
		continue;
	pSrc = szLineBuff;
	pEntryFld = pEntryOfsFld = pValFld = NULL;
	FldPsn = 1;
	if(FldPsn == EntryFld)
		pEntryFld = pSrc;
	if(FldPsn == EntryOfsFld)
		pEntryOfsFld = pSrc;
	if(FldPsn == ValFld)
		pValFld = pSrc;
	while((chr = *pSrc++)!= '\0')
		{
		if(chr == ',')
			{
			FldPsn++;
			if(FldPsn == EntryFld)
				pEntryFld = pSrc;
			if(FldPsn == EntryOfsFld)
				pEntryOfsFld = pSrc;
			if(FldPsn == ValFld)
				pValFld = pSrc;
			}
		}
	if(pEntryFld == NULL)
		{
		AddErrMsg("CCSV2BED::process","Unable to locate the entry field at %d in file %s",EntryFld,pszCSV);
		fclose(pCSVStream);
		return(eBSFerrLocField);
		}
	if(pEntryOfsFld == NULL)
		{
		AddErrMsg("CCSV2BED::process","Unable to locate the entry offset field at %d in file %s",EntryOfsFld,pszCSV);
		fclose(pCSVStream);
		return(eBSFerrLocField);
		}
	if(pValFld == NULL)
		{
		AddErrMsg("CCSV2BED::process","Unable to locate the value field at %d in file %s",ValFld,pszCSV);
		fclose(pCSVStream);
		return(eBSFerrLocField);
		}

	if(sscanf(pEntryFld,"%d,",&Entry)!=1)
		{
		AddErrMsg("CCSV2BED::process","Unable to parse the entry field at %d as a integer in file %s",EntryFld,pszCSV);
		fclose(pCSVStream);
		return(eBSFerrLocField);
		}

	if(Entry != EntryID)
		continue;

	if(sscanf(pEntryOfsFld,"%d",&CurEntryOfs)!=1)
		{
		AddErrMsg("CCSV2BED::process","Unable to parse the entry offset field at %d as a integer in file %s",EntryOfsFld,pszCSV);
		fclose(pCSVStream);
		return(eBSFerrLocField);
		}

	// filter out entries which are outside the range specified 
	if(CurEntryOfs < StartOfs || (EndOfs && (CurEntryOfs > EndOfs)))
		continue;
	if(sscanf(pValFld,"%lf",&StepValue)!=1)
		{
		AddErrMsg("CCSV2BED::process","Unable to parse the value field at %d as a float in file %s",EntryFld,pszCSV);
		fclose(pCSVStream);
		return(eBSFerrLocField);
		}

	AddStepValue(CurEntryOfs + Remap,(int)(StepValue * 10000));
	}
fclose(pCSVStream);

if(!m_NumSteps)
	{
	AddErrMsg("CCSV2BED::process","No steps located for entry %d between %d and %d in %s",EntryID,StartOfs,EndOfs,pszCSV);
	return(eBSFerrLocField);
	}

// CSV may not have had entry offsets in ascending order so do a quick sort on them
qsort(m_pStepValues,m_NumSteps,sizeof(tsStepValues),CompareStepPsns);

// start data reduction and scaling ready for BED output
return(Output2BEDorWIG(pszBEDorWIG,pszTrackName,pszDescr,pszChrom,pValName,NumClasses,Windowsize,bAppend,bWIG));
}

int
CCSV2BED::AddStepValue(int StepPsn,long StepValue)
{
tsStepValues *pNewAlloc;
int NumSteps2Alloc;

	// need to grow m_pStepValues allocation?
if(m_pStepValues == NULL || m_NumSteps >= m_NumStepsAllocd)
	{
	NumSteps2Alloc = m_NumStepsAllocd + 10000000;
	pNewAlloc = (tsStepValues *) new char [sizeof(tsStepValues) * NumSteps2Alloc];
	if(pNewAlloc == NULL)
		return(eBSFerrMem);
	if(m_pStepValues != NULL)
		{
		if(m_NumSteps)
			memmove(pNewAlloc,m_pStepValues,m_NumSteps * sizeof(tsStepValues));
		delete m_pStepValues;
		}
	else
		m_NumSteps = 0;
	m_pStepValues = pNewAlloc;
	m_NumStepsAllocd = NumSteps2Alloc;
	}
m_pStepValues[m_NumSteps].StepPsn = StepPsn;
m_pStepValues[m_NumSteps++].StepValue = StepValue;
return(m_NumSteps);
}

// Yk = mean for current sample
// Yj = previous mean
// N  = window size
// Xk = current sample
// Xj = previous sample
// Yk = Yj +  1/N * (Xk - Xj)
// 
bool
CCSV2BED::SmoothData(int WindowSize)
{
if(WindowSize < 2 || m_pStepValues == NULL || m_NumSteps < WindowSize)
	return(false);
if(WindowSize > 100)
	WindowSize = 100;

tsStepValues *pYk;
int CurStep;
int N,Xj,Yk,Xk;

pYk = m_pStepValues;
Xj = pYk->StepValue;
Yk = pYk->StepValue;
N = 1;
for(CurStep = 0; CurStep < m_NumSteps; CurStep++,pYk++, N >= WindowSize ? N : N++)
	{
	Xk = Xj;
	Xj = pYk->StepValue;
	Yk += (Xj - Xk)/N;
	pYk->StepValue = Yk;
	}
return(true);
}

// Output2BEDorWIG
// Output steps + values to BED or WIGGLE file
// BED format --
// <chrom> <chromStart (0..n)> <chromEnd (chromStart+Len)> <name> <score (0..999)>
// chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or contig (e.g. ctgY1). 
// chromStart - The starting position of the feature in the chromosome or contig. The first base in a chromosome is numbered 0. 
// chromEnd - The ending position of the feature in the chromosome or contig. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
// name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode. 
// score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). 
//
// WIG format --
// As BED format:
// <chrom>  <chromStart (1..n)>  <chromEnd (chromStart+Len)>  <dataValue>


int
CCSV2BED::Output2BEDorWIG(char *pszBED,		// BED or WIG file to create/append
					 char *pszTrackName,	// track name="???"
					 char *pszDescr,		// description="???"
					 char *pszChrom,		// name of chromosome
					 char *pValName,		// name to output as BED or WIG
					 int NumClasses,		// number of classes to bin mean values into
					int Windowsize,			// window size to use for data mean smoothing (1 if no smoothing)
					 bool bAppend,			// true to append on to existing file
					 bool bWIG)				// true if to output as wiggle file, default is as BED
{
char szLineBuff[200];
int LineLen;
int hBED;
int CurStep;
int MinStepValue;
int MaxStepValue;
int StepValueRange;
int ClassRange;
int CurStartStep;
tsStepValues *pCurStep;
char szVal[80];
int openopts;
int MeanWinLen;
int CurMeanScore;
int MeanClassSize;
int CurMean;
int NewMean;
int MinMean;
int MaxMean;
unsigned int ExpStep;

if(m_pStepValues == NULL || m_NumSteps == 0)
	return(eBSFerrParams);

ClassRange = 1000 / NumClasses;			// BED max score is 1000, so determine range of each class

if(Windowsize > 1)
	SmoothData(Windowsize);				// apply data smoothing (mean over sliding window of Windowsize)

// determine min/max values, and hence StepValueRange
MaxStepValue = 0x80000000;
MinStepValue = 0x7fffffff;
pCurStep = m_pStepValues;
for(CurStep = 0; CurStep < m_NumSteps; CurStep++,pCurStep++)
	{
	if(pCurStep->StepValue > MaxStepValue)
		MaxStepValue = pCurStep->StepValue;
	if(pCurStep->StepValue < MinStepValue)
		MinStepValue = pCurStep->StepValue;
	}
StepValueRange = MaxStepValue - MinStepValue;
MeanClassSize = (StepValueRange + NumClasses) / NumClasses;

// create or open existing output file taking into account if the file is to be appended or truncated
#ifdef _WIN32
openopts = _O_RDWR | _O_BINARY | _O_SEQUENTIAL | _O_CREAT;
openopts |= bAppend ? 0 : _O_TRUNC;
hBED = open(pszBED,openopts, _S_IREAD | _S_IWRITE );
#else
openopts = O_RDWR | O_CREAT;
openopts |= bAppend ? 0 : O_TRUNC;
hBED = open(pszBED,openopts, S_IREAD | S_IWRITE );
#endif
if(hBED == -1)
	{
	AddErrMsg("CCSV2BED::Output2BEDorWIG","Unable to open %s - %s",pszBED,strerror(errno));
	return(eBSFerrCreateFile);
	}
if(bAppend)
	_lseeki64(hBED,0,SEEK_END);

if(bWIG)
	LineLen = sprintf(szLineBuff,"\ntrack type=wiggle_0 autoscale=on windowingFunction=mean smoothingWindow=5 name=%s description=\"%s\" useScore=1\n",pszTrackName,pszDescr);
else
	LineLen = sprintf(szLineBuff,"\ntrack name=%s description=\"%s\" useScore=1\n",pszTrackName,pszDescr);

if(write(hBED,szLineBuff,LineLen)!=LineLen)
	{
	AddErrMsg("CCSV2BED::Output2BEDorWIG","Write to %s - %s",pszBED,strerror(errno));
	close(hBED);
	return(eBSFerrFileAccess);
	}

// iterate over all steps
// if the step value is outside the current class bin range, or if there is a gap in steps, then output the current
// class bin mean value and start a new bin
pCurStep = m_pStepValues;
CurStartStep = pCurStep->StepPsn;
ExpStep = CurStartStep + 1;
CurMean = pCurStep->StepValue;
MinMean = MinStepValue + ((CurMean - MinStepValue)/MeanClassSize)*MeanClassSize;
MaxMean = MinMean + MeanClassSize;
MeanWinLen = 1;
pCurStep++;
for(CurStep = 2; CurStep < m_NumSteps; CurStep++,pCurStep++,MeanWinLen++,ExpStep++)
	{
	NewMean = pCurStep->StepValue;

	if(ExpStep == pCurStep->StepPsn &&					// steps (psn) has to be continuous  
		(NewMean >= MinMean && NewMean <= MaxMean))			// and within same mean class range to continue
		{
		CurMean += NewMean;
		continue;
		}

	CurMean /= MeanWinLen;
	if(bWIG) // is WIGGLE format
		// <chrom>  <chromStart (1..n)>  <chromEnd (chromStart+Len)>  <dataValue>
		LineLen = sprintf(szLineBuff,"%s %d %d %s\n",pszChrom,CurStartStep+1,ExpStep+1,Fmt2FixedDec(CurMean,szVal));
	else // is BED format
		{
		// <chrom> <chromStart (0..n)> <chromEnd (chromStart + Len)> <name> <score (0..999) >
		CurMeanScore = ((CurMean - MinStepValue) * 1000)/StepValueRange;
		LineLen = sprintf(szLineBuff,"%s %d %d %s %d\n",pszChrom,CurStartStep,ExpStep,
					pValName[0] != '\0' ? pValName : Fmt2FixedDec(CurMean,szVal), // use mean value as the name
					CurMeanScore);
		}

	if(write(hBED,szLineBuff,LineLen)!=LineLen)
		{
		AddErrMsg("CCSV2BED::Output2BEDorWIG","Write to %s - %s",pszBED,strerror(errno));
		close(hBED);
		return(eBSFerrFileAccess);
		}
	CurMean = NewMean;
	ExpStep = CurStartStep = pCurStep->StepPsn;
	MeanWinLen = 0;
	MinMean = MinStepValue + ((CurMean - MinStepValue)/MeanClassSize)*MeanClassSize;
	MaxMean = MinMean + MeanClassSize;
	}

CurMean /= MeanWinLen;
if(bWIG)
		LineLen = sprintf(szLineBuff,"%s %d %d %s\n",pszChrom,CurStartStep+1,ExpStep+1,Fmt2FixedDec(CurMean,szVal));
else
	{
	CurMeanScore = ((CurMean - MinStepValue) * 1000)/StepValueRange;
	LineLen = sprintf(szLineBuff,"%s %d %d %s %d\n",pszChrom,CurStartStep,ExpStep,
			pValName[0] != '\0' ? pValName : Fmt2FixedDec(CurMean,szVal),  // use mean value as the name
			CurMeanScore);
	}

if(write(hBED,szLineBuff,LineLen)!=LineLen)
	{
	AddErrMsg("CCSV2BED::Output2BEDorWIG","Write to %s - %s",pszBED,strerror(errno));
	close(hBED);
	return(eBSFerrFileAccess);
	}
close(hBED);
return(eBSFSuccess);
}


int CCSV2BED::CompareStepPsns( const void *arg1, const void *arg2 )
{
tsStepValues *pEl1 = (tsStepValues *)arg1;
tsStepValues *pEl2 = (tsStepValues *)arg2;
if(pEl1->StepPsn == pEl2->StepPsn)
	return(0);
return(pEl1->StepPsn < pEl2->StepPsn ? -1 : 1);
}

char *
CCSV2BED::Fmt2FixedDec(int Value,char *pszRet)
{
int Intg = Value/10000;
int Remainder = abs(Value)%10000;
sprintf(pszRet,"%d.%04d",Intg,Remainder);
return(pszRet);
}

