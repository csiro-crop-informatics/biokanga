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


char CSeqTrans::MapBase2Ascii(etSeqBase Base,		// base to map 
							char BaseN ,			// what to map eBaseN onto
							char BaseUndef,			// what to map eBaseUndef onto
							char BaseInDel,			// what to map eBaseInDel onto
							bool RptMskUpperCase)
{
if(RptMskUpperCase) 
	{
	switch(Base  & ~cMarkMskFlg) {
			case eBaseA: return('a');
			case eBaseA | cRptMskFlg: return('A');

			case eBaseC: return('c');
			case eBaseC | cRptMskFlg: return('C');

			case eBaseG: return('g');
			case eBaseG | cRptMskFlg: return('G');

			case eBaseT: return('t');
			case eBaseT | cRptMskFlg: return('T');

			case eBaseN: return(tolower(BaseN));
			case eBaseN | cRptMskFlg: return(toupper(BaseN));

			case eBaseInDel: return(BaseInDel);
			case eBaseUndef: return(BaseUndef);
			default: return('?');
			}
	}
else
	{
	switch(Base  & ~cMarkMskFlg) {
			case eBaseA: return('A');
			case eBaseA | cRptMskFlg: return('a');

			case eBaseC: return('C');
			case eBaseC | cRptMskFlg: return('c');

			case eBaseG: return('G');
			case eBaseG | cRptMskFlg: return('g');

			case eBaseT: return('T');
			case eBaseT | cRptMskFlg: return('t');

			case eBaseN: return(toupper(BaseN));
			case eBaseN | cRptMskFlg: return(tolower(BaseN));

			case eBaseInDel: return(BaseInDel);
			case eBaseUndef: return(BaseUndef);
			default: return('?');
			}
	}
}

char *
CSeqTrans::MapSeq2Ascii(etSeqBase *pBases,			// bases to map to ascii
							unsigned int SeqLen,	// number of bases to map
							char *pAscii,			// where to return ascii - NOTE: will have terminating '\0' appended
							char BaseN,				// what to map eBaseN onto
							char BaseUndef,			// what to map eBaseUndef onto
							char BaseInDel,			// what to map eBaseInDel onto
							bool RptMskUpperCase)   // if true then uppercase represents softmasked repeats
{
static char szSeq2AsciiBuff[cMaxAutoSeq2Ascii+1];
char *pszBuffer;
if(pAscii != NULL)
	pszBuffer = pAscii;
else
	{
	pAscii = szSeq2AsciiBuff;
	pszBuffer = szSeq2AsciiBuff;
	if(SeqLen >= sizeof(szSeq2AsciiBuff))
		SeqLen = sizeof(szSeq2AsciiBuff)-1;
	}
while(SeqLen--)
	{
	if(RptMskUpperCase) {

		switch(*pBases++ & ~cMarkMskFlg) {
			case eBaseA: *pszBuffer++ = 'a'; continue;
			case eBaseA | cRptMskFlg: *pszBuffer++ = 'A'; continue;

			case eBaseC: *pszBuffer++ = 'c'; continue;
			case eBaseC | cRptMskFlg: *pszBuffer++ = 'C'; continue;

			case eBaseG: *pszBuffer++ = 'g'; continue;
			case eBaseG | cRptMskFlg: *pszBuffer++ = 'G'; continue;

			case eBaseT: *pszBuffer++ = 't'; continue;
			case eBaseT | cRptMskFlg: *pszBuffer++ = 'T'; continue;

			case eBaseN: *pszBuffer++ = tolower(BaseN); continue;
			case eBaseN | cRptMskFlg: *pszBuffer++ = toupper(BaseN); continue;

			case eBaseInDel: *pszBuffer++ = BaseInDel; continue;
			case eBaseUndef: *pszBuffer++ = BaseUndef; continue;
			default: 
				*pszBuffer++ = '?'; 
				continue;
			}
		}
	else
		{
		switch(*pBases++ & ~cMarkMskFlg) {
			case eBaseA: *pszBuffer++ = 'A'; continue;
			case eBaseA | cRptMskFlg: *pszBuffer++ = 'a'; continue;

			case eBaseC: *pszBuffer++ = 'C'; continue;
			case eBaseC | cRptMskFlg: *pszBuffer++ = 'c'; continue;

			case eBaseG: *pszBuffer++ = 'G'; continue;
			case eBaseG | cRptMskFlg: *pszBuffer++ = 'g'; continue;

			case eBaseT: *pszBuffer++ = 'T'; continue;
			case eBaseT | cRptMskFlg: *pszBuffer++ = 't'; continue;

			case eBaseN: *pszBuffer++ = toupper(BaseN); continue;
			case eBaseN | cRptMskFlg: *pszBuffer++ = tolower(BaseN); continue;

			case eBaseInDel: *pszBuffer++ = BaseInDel; continue;
			case eBaseUndef: *pszBuffer++ = BaseUndef; continue;
			default: 
				*pszBuffer++ = '?'; 
				continue;
			}
		}
	break;
	}
*pszBuffer = '\0';
return(pAscii);
}

// MapSeq2UCAscii
// Returns ascii sequence with all bases as uppercase
char *
CSeqTrans::MapSeq2UCAscii(etSeqBase *pBases,			// bases to map to ascii
							unsigned int SeqLen,	// number of bases to map
							char *pAscii,			// where to return ascii - NOTE: will have terminating '\0' appended
							char BaseN,				// what to map eBaseN onto
							char BaseUndef,			// what to map eBaseUndef onto
							char BaseInDel)			// what to map eBaseInDel onto
{
static char szSeq2AsciiBuff[cMaxAutoSeq2Ascii+1];
char *pszBuffer;
if(pAscii != NULL)
	pszBuffer = pAscii;
else
	{
	pAscii = szSeq2AsciiBuff;
	pszBuffer = szSeq2AsciiBuff;
	if(SeqLen >= sizeof(szSeq2AsciiBuff))
		SeqLen = sizeof(szSeq2AsciiBuff)-1;
	}
while(SeqLen--)
	{
	switch(*pBases++ & 0x07) {
			case eBaseA: *pszBuffer++ = 'A'; continue;
			case eBaseC: *pszBuffer++ = 'C'; continue;
			case eBaseG: *pszBuffer++ = 'G'; continue;
			case eBaseT: *pszBuffer++ = 'T'; continue;
			case eBaseN: *pszBuffer++ = BaseN; continue;
			case eBaseInDel: *pszBuffer++ = BaseInDel; continue;
			case eBaseUndef: *pszBuffer++ = BaseUndef; continue;
			default: 
				*pszBuffer++ = '?'; continue;
		}
	break;
	}
*pszBuffer = '\0';
return(pAscii);
}

// MapSeq2LCAscii
// Returns ascii sequence with all bases as lowercase
char *
CSeqTrans::MapSeq2LCAscii(etSeqBase *pBases,			// bases to map to ascii
							unsigned int SeqLen,	// number of bases to map
							char *pAscii,			// where to return ascii - NOTE: will have terminating '\0' appended
							char BaseN,				// what to map eBaseN onto
							char BaseUndef,			// what to map eBaseUndef onto
							char BaseInDel)			// what to map eBaseInDel onto
{
static char szSeq2AsciiBuff[cMaxAutoSeq2Ascii+1];
char *pszBuffer;
if(pAscii != NULL)
	pszBuffer = pAscii;
else
	{
	pAscii = szSeq2AsciiBuff;
	pszBuffer = szSeq2AsciiBuff;
	if(SeqLen >= sizeof(szSeq2AsciiBuff))
		SeqLen = sizeof(szSeq2AsciiBuff)-1;
	}
while(SeqLen--)
	{
	switch(*pBases++ & 0x07) {
			case eBaseA: *pszBuffer++ = 'a'; continue;
			case eBaseC: *pszBuffer++ = 'c'; continue;
			case eBaseG: *pszBuffer++ = 'g'; continue;
			case eBaseT: *pszBuffer++ = 't'; continue;
			case eBaseN: *pszBuffer++ = BaseN; continue;
			case eBaseInDel: *pszBuffer++ = BaseInDel; continue;
			case eBaseUndef: *pszBuffer++ = BaseUndef; continue;
			default: 
				*pszBuffer++ = '?'; continue;
		}
	break;
	}
*pszBuffer = '\0';
return(pAscii);
}



char *
CSeqTrans::MapPackedSeq2Ascii(int StartNibble,	// which nibble to start from - 0 if lower or first nibble
							etSeqBase *pBases,	// packed bases (bases in lower and upper nibbles) to map to ascii
							unsigned int SeqLen,	// total number of bases to map
							char *pAscii,			// where to return ascii - NOTE: will have terminating '\0' appended
							char BaseN,				// what to map eBaseN onto
							char BaseUndef,			// what to map eBaseUndef onto
							char BaseInDel,			// what to map eBaseInDel onto
							bool RptMskUpperCase)   // if true then uppercase represents softmasked repeats
{
etSeqBase SeqBase;
static char szSeq2AsciiBuff[cMaxAutoSeq2Ascii+1];
char *pszBuffer;

if(pAscii != NULL)
	pszBuffer = pAscii;
else
	{
	pAscii = szSeq2AsciiBuff;
	pszBuffer = szSeq2AsciiBuff;
	if(SeqLen >= sizeof(szSeq2AsciiBuff))
		SeqLen = sizeof(szSeq2AsciiBuff)-1;
	}

while(SeqLen--)
	{
	if(!(StartNibble & 0x01))
		SeqBase = *pBases & 0x0f;
	else
		SeqBase = (*pBases++ >> 4) & 0x0f;
	StartNibble += 1;
	if(RptMskUpperCase) {
		switch(SeqBase) {
			case eBaseA: *pszBuffer++ = 'a'; continue;
			case eBaseA | cRptMskFlg: *pszBuffer++ = 'A'; continue;

			case eBaseC: *pszBuffer++ = 'c'; continue;
			case eBaseC | cRptMskFlg: *pszBuffer++ = 'C'; continue;

			case eBaseG: *pszBuffer++ = 'g'; continue;
			case eBaseG | cRptMskFlg: *pszBuffer++ = 'G'; continue;

			case eBaseT: *pszBuffer++ = 't'; continue;
			case eBaseT | cRptMskFlg: *pszBuffer++ = 'T'; continue;

			case eBaseN: *pszBuffer++ = tolower(BaseN); continue;
			case eBaseN | cRptMskFlg: *pszBuffer++ = toupper(BaseN); continue;

			case eBaseInDel: *pszBuffer++ = BaseInDel; continue;
			case eBaseUndef: *pszBuffer++ = BaseUndef; continue;
			default: *pszBuffer++ = '?'; break;
			}
		}
	else
		{
		switch(SeqBase) {
			case eBaseA: *pszBuffer++ = 'A'; continue;
			case eBaseA | cRptMskFlg: *pszBuffer++ = 'a'; continue;

			case eBaseC: *pszBuffer++ = 'C'; continue;
			case eBaseC | cRptMskFlg: *pszBuffer++ = 'c'; continue;

			case eBaseG: *pszBuffer++ = 'G'; continue;
			case eBaseG | cRptMskFlg: *pszBuffer++ = 'g'; continue;

			case eBaseT: *pszBuffer++ = 'T'; continue;
			case eBaseT | cRptMskFlg: *pszBuffer++ = 't'; continue;

			case eBaseN: *pszBuffer++ = toupper(BaseN); continue;
			case eBaseN | cRptMskFlg: *pszBuffer++ = tolower(BaseN); continue;

			case eBaseInDel: *pszBuffer++ = BaseInDel; continue;
			case eBaseUndef: *pszBuffer++ = BaseUndef; continue;
			default: *pszBuffer++ = '?'; break;
			}
		}
	break;
	}
*pszBuffer = '\0';
return(pAscii);
}


etSeqBase *
CSeqTrans::MapAscii2Sense(char *pAscii,
							unsigned int AsciiLen,
							etSeqBase *pSense,
							bool RptMskUpperCase)   // if true then uppercase represents softmasked repeats
{
static unsigned char sAscii2SeqBuff[cMaxAutoSeq2Ascii+1];
unsigned char *pSeq;

if(AsciiLen < 0 || pAscii == NULL || *pAscii == '\0')
	return(NULL);
if(AsciiLen == 0)
	AsciiLen = (unsigned int)strlen(pAscii);
if(pSense != NULL)
	pSeq = pSense;
else
	{
	pSeq = sAscii2SeqBuff;
	pSense = pSeq;
	if(AsciiLen >= sizeof(sAscii2SeqBuff))
		AsciiLen = sizeof(sAscii2SeqBuff)-1;
	}

while(AsciiLen--)
	{
	if(RptMskUpperCase)
		{
		switch(*pAscii++) {
			case ' ': continue;					// simply slough any contained spaces, tabs, or the comma being used as a separator 
			case '\t': continue;
			case ',': continue;
			case 'a': *pSeq++ = eBaseA; continue;
			case 'A': *pSeq++ = eBaseA | cRptMskFlg;	continue;
			case 'c': *pSeq++ = eBaseC;	continue;
			case 'C': *pSeq++ = eBaseC | cRptMskFlg;	continue;
			case 'g':*pSeq++ = eBaseG; continue;
			case 'G':*pSeq++= eBaseG | cRptMskFlg; continue;
			case 't': case 'u': *pSeq++ = eBaseT;	continue;
			case 'T': case 'U':*pSeq++ = eBaseT | cRptMskFlg;	continue;
			case 'n': *pSeq++ = eBaseN; continue;
			case 'N': *pSeq++ = eBaseN | cRptMskFlg;	continue;
			case '-':*pSeq++ = eBaseInDel; continue;
			case '\0': *pSeq = eBaseEOS; return(pSense);
			default: *pSeq++ = eBaseUndef; continue;
			}
		}
	else
		{
		switch(*pAscii++) {
			case ' ': continue;					// simply slough any contained spaces, tabs, or the comma being used as a separator 
			case '\t': continue;
			case ',': continue;
			case 'A': *pSeq++ = eBaseA; continue;
			case 'a': *pSeq++ = eBaseA | cRptMskFlg;	continue;
			case 'C': *pSeq++ = eBaseC;	continue;
			case 'c': *pSeq++ = eBaseC | cRptMskFlg;	continue;
			case 'G':*pSeq++ = eBaseG; continue;
			case 'g':*pSeq++= eBaseG | cRptMskFlg; continue;
			case 'T': case 'U': *pSeq++ = eBaseT;	continue;
			case 't': case 'u':*pSeq++ = eBaseT | cRptMskFlg;	continue;
			case 'n': *pSeq++ = eBaseN | cRptMskFlg;	continue;
			case 'N': *pSeq++ = eBaseN;	continue;
			case '-':*pSeq++ = eBaseInDel; continue;
			case '\0': *pSeq = eBaseEOS; return(pSense);
			default: *pSeq++ = eBaseUndef; continue;
			}
		}
	}
return(pSense);
}


etSeqBase * 
CSeqTrans::MapAscii2Antisense(char *pAscii,
							unsigned int AsciiLen,
							etSeqBase *pAntisense,
							bool RptMskUpperCase)   // if true then uppercase represents softmasked repeats
{
static unsigned char sAscii2SeqBuff[cMaxAutoSeq2Ascii+1];
unsigned char *pSeq;

if(AsciiLen < 0 || pAscii == NULL || *pAscii == '\0')
	return(NULL);
if(AsciiLen == 0)
	AsciiLen = (unsigned int)strlen(pAscii);
if(pAntisense != NULL)
	pSeq = pAntisense;
else
	{
	pSeq = sAscii2SeqBuff;
	pAntisense = pSeq;
	if(AsciiLen >= sizeof(sAscii2SeqBuff))
		AsciiLen = sizeof(sAscii2SeqBuff)-1;
	}

while(AsciiLen--)
	{
	if(RptMskUpperCase)
		{
		switch(*pAscii++) {
			case ' ': continue;					// simply slough any contained spaces, tabs, or the comma being used as a separator 
			case '\t': continue;
			case ',': continue;
			case 'a': *pSeq++ = eBaseT; continue;
			case 'A': *pSeq++ = eBaseT | cRptMskFlg;	continue;
			case 'c': *pSeq++ = eBaseG;	continue;
			case 'C': *pSeq++ = eBaseG | cRptMskFlg;	continue;
			case 'g':*pSeq++ = eBaseC; continue;
			case 'G':*pSeq++= eBaseC | cRptMskFlg; continue;
			case 't': case 'u': *pSeq++ = eBaseA;	continue;
			case 'T': case 'U':*pSeq++ = eBaseA | cRptMskFlg;	continue;
			case '-':*pSeq++ = eBaseInDel; continue;
			case 'n':*pSeq++ = eBaseN; continue;
			case 'N':*pSeq++ = eBaseN  | cRptMskFlg; continue;
			case '\0': *pSeq = eBaseEOS; return(pAntisense);
			default: *pSeq++ = eBaseUndef; continue;
			}
		}
	else
		{
		switch(*pAscii++) {
			case ' ': continue;					// simply slough any contained spaces, tabs, or the comma being used as a separator 
			case '\t': continue;
			case ',': continue;
			case 'A': *pSeq++ = eBaseT; continue;
			case 'a': *pSeq++ = eBaseT | cRptMskFlg;	continue;
			case 'C': *pSeq++ = eBaseG;	continue;
			case 'c': *pSeq++ = eBaseG | cRptMskFlg;	continue;
			case 'G':*pSeq++ = eBaseC; continue;
			case 'g':*pSeq++= eBaseC | cRptMskFlg; continue;
			case 'T': case 'u': *pSeq++ = eBaseA;	continue;
			case 't': case 'U':*pSeq++ = eBaseA | cRptMskFlg;	continue;
			case '-':*pSeq++ = eBaseInDel; continue;
			case 'N':*pSeq++ = eBaseN; continue;
			case 'n':*pSeq++ = eBaseN  | cRptMskFlg; continue;
			case '\0': *pSeq = eBaseEOS; return(pAntisense);
			default: *pSeq++ = eBaseUndef; continue;
			}
		}
	}
return(pAntisense);
}


void 
CSeqTrans::ComplementStrand(unsigned int SeqLen,
							etSeqBase *pSeq)
{
etSeqBase Base;
etSeqBase RptMskFlg;
if(SeqLen < 1 || pSeq == NULL)
	return;
while(SeqLen--)
	{
	Base = *pSeq;
	RptMskFlg = Base & (cRptMskFlg | cMarkMskFlg);
	Base &= ~(cRptMskFlg | cMarkMskFlg);
	switch(Base) {
		case eBaseA: *pSeq++ = eBaseT | RptMskFlg; continue;
		case eBaseC: *pSeq++ = eBaseG | RptMskFlg; continue;
		case eBaseG: *pSeq++ = eBaseC | RptMskFlg; continue;
		case eBaseT: *pSeq++ = eBaseA | RptMskFlg; continue;
		case eBaseN: case eBaseInDel: case eBaseUndef: pSeq++; continue;
		default:
			return;
		}
	}
}


// ReverseSeq
// Inplace sequence reversal agtc becomes ctga
void
CSeqTrans::ReverseSeq(unsigned int SeqLen,etSeqBase *pSeq)
{
etSeqBase *pExch;
etSeqBase Tmp;
if(SeqLen < 2 || pSeq == NULL)
	return;
pExch = pSeq + SeqLen-1;
SeqLen >>= 1;

while(SeqLen--)
	{
	Tmp = *pExch;
	*pExch-- = *pSeq;
	*pSeq++ = Tmp;
	}
}

// ReverseComplement
// Inplace sequence reverse complement agtc becomes gact
void
CSeqTrans::ReverseComplement(unsigned int SeqLen,etSeqBase *pSeq)
{
if(SeqLen < 1 || pSeq == NULL)
	return;
ComplementStrand(SeqLen,pSeq);
ReverseSeq(SeqLen,pSeq);
}

// RemoveMasking
// Removes all softmasking from specified sequence
void 
CSeqTrans::RemoveMasking(etSeqBase *pBases,int SeqLen)
{
if(pBases == NULL || SeqLen <= 0)
	return;
while(SeqLen--)
	*pBases++ &= ~(cRptMskFlg | cMarkMskFlg);
}

