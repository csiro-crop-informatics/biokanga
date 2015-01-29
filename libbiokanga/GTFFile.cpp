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
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

typedef struct TAG_sFeatTypes {
	etGGTFFeatType FeatType;
	const char *pszType;
	const char *pszDescr;
} tsFeatType;

static tsFeatType m_FeatTypes[] = {
	{eGGTFCDS,"CDS","coding sequence"},
	{eGGTFStartCodon,"start_codon","start codon"},
	{eGGTFStopCodon,"stop_codon","stop codon"},
	{eGGTF5UTR,"5UTR","5' UTR"},
	{eGGTF3UTR,"3UTR","3' UTR"},
	{eGGTFinter,"inter","intergenic region"},
	{eGGTFinter_CNS,"inter_CNS","intergenic conserved noncoding sequence region"},
	{eGGTFintron_CNS,"intron_CNS","conserved noncoding sequence region within an intron"},
	{eGGTFexon,"exon","exon"}
};
int m_NumFeatTypes = (sizeof(m_FeatTypes)/sizeof(tsFeatType));

CGTFFile::CGTFFile(void)
{
m_pGTFStream = NULL;
Reset();
}

CGTFFile::~CGTFFile(void)
{
Reset();
}

int 
CGTFFile::Open(char *pszFileName)		// file to open
{
Reset();

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((m_pGTFStream = fopen(pszFileName,"r"))==NULL)
	{
	AddErrMsg("CGTFFile::Open","Unable to open GTF format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
strcpy(m_szInFile,pszFileName);
memset(&m_Fields,0,sizeof(m_Fields));
m_bEOF = false;
return(eBSFSuccess);
}

int 
CGTFFile::Close(void)
{
Reset();
return(eBSFSuccess);
}

void
CGTFFile::Reset(void)
{
if(m_pGTFStream != NULL)
	{
	fclose(m_pGTFStream);
	m_pGTFStream = NULL;
	}

m_szInFile[0]='\0';

memset(&m_Fields,0,sizeof(m_Fields));
m_CurLineNum = 0;		// current line number
m_CurNumFields = 0;		// number of fields parsed from current line
m_CurMaxFieldLen = 0;	// max len of any field in current line
m_bEOF = false;			// TRUE after last record parsed
}

char *
CGTFFile::TrimWhitespace(char *pTxt)
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

char *
CGTFFile::GetNoteValue(char *pszAttribute)
{
int TagLen;
char *pszAttribTag;
char *pszAttribVal;
static char szValue[cGTFMaxValLen];
if(pszAttribute == NULL || m_pGTFStream == NULL || m_bEOF || m_Fields.AttribOfs == 0)
	return(NULL);
szValue[0] = '\0';
TagLen = (int)strlen(pszAttribute);
pszAttribTag = &m_Fields.szRawLine[m_Fields.AttribOfs];
while(*pszAttribTag != '\0')
	{
	if(*pszAttribTag == '\t' || *pszAttribTag == ' ' || *pszAttribTag == ';')
		{
		pszAttribTag += 1;
		continue;
		}
	if(!strnicmp(pszAttribute,pszAttribTag,TagLen))
		{
		pszAttribTag += TagLen + 1;
		pszAttribVal = szValue;
		TagLen = cGTFMaxValLen-1;
		while((*pszAttribVal++ = *pszAttribTag++) && TagLen--)
			if(*pszAttribTag == ';')
				break;
		*pszAttribVal = '\0';
		return(szValue);
		}
	// check against next attribute
	while(*pszAttribTag != '\0' && *pszAttribTag != ';')
		pszAttribTag += 1;
	}
return(NULL);
}

// Iterate records associated with specified type
// "start_codon","stop_codon","exon","CDS", etc
// This requires a check on both the feature field and also the Note attribute as well as
// tracking on the Parent attribute
//
int					// 1 if record accepted, 0 if reached EOF, < 0 if error
CGTFFile::NextRecordOfType(etGGTFFeatType GGTFFeatType)
{
int Rslt;
int Idx;
tsFeatType *pFeatType;
static char szCurID[cMaxGeneNameLen] = "\0";
char *pszNoteVal = NULL;


if(GGTFFeatType < eGGTFany || GGTFFeatType >= eGGTFplaceholder)
	{
	AddErrMsg("CGTFFile::NextRecordOfType","Requested GTF type '%d' not supported",GGTFFeatType);
	return(eBSFerrFeature);
	}

if(GGTFFeatType != eGGTFany)
	{
	pFeatType = m_FeatTypes;
	for(Idx = 0; Idx < m_NumFeatTypes; Idx++,pFeatType++)
		if(pFeatType->FeatType == GGTFFeatType)
			break;
	}

while((Rslt = NextLine()) > 0)
	{
	if(GGTFFeatType == eGGTFany)
		return(1);
	if(!stricmp(m_Fields.szFeature,pFeatType->pszType))
		return(1);
	}
return(Rslt);
}

int							// 1 if record processed, 0 if reached EOF, < 0 if error
CGTFFile::NextLine(void)			// move to next line in GTF file, skips blank lines and comment lines starting with comment char
{
return(ParseLine());
}

int 
CGTFFile::ParseLine(void)	// read next line and parse
{
int Cnt;
char *pTxt;
int AttribListStart;
char szScore[50];
char cFrame;

if(m_pGTFStream == NULL)
	{
	AddErrMsg("CGTFFile::ParseLine","GTF input stream closed");
	return(eBSFerrFileClosed);
	}
if(m_bEOF)				// last record already parsed?
	return(0);
m_Fields.AttribOfs = 0;
m_Fields.End = 0;
m_Fields.bDfltStrand = true;
m_Fields.Frame = 0;
m_Fields.bDfltScore = true;
m_Fields.Score = 0.0f;
m_Fields.Start = 0;
m_Fields.ScoreOfs = 0;
m_Fields.bDfltStrand = true;
m_Fields.Strand = '?';
m_Fields.szFeature[0] = '\0';
m_Fields.szRawLine[0] = '\0';
m_Fields.szSeqName[0] = '\0';
m_Fields.szSeqSource[0] = '\0';
AttribListStart = 0;

while(fgets(m_Fields.szRawLine,sizeof(m_Fields.szRawLine)-1,m_pGTFStream)!=NULL)
	{
	m_CurLineNum += 1;
	pTxt = TrimWhitespace(m_Fields.szRawLine);
	m_CurLineLen = (int)strlen(pTxt);
	if(*pTxt=='\0' || (*pTxt=='#' && (pTxt[1] != '#' || pTxt[1] == '\0')))	// simply slough lines which were just whitespace or start with '#'
		continue;
	if(pTxt[0] == '#' && pTxt[1] == '#')
		{
		return(1);
		}
	// parse out the mandatory fields
	Cnt = sscanf(pTxt,"%34s\t%34s\t%34s\t%d\t%d\t%n%25s\t%c\t%c%n",
			m_Fields.szSeqName,m_Fields.szSeqSource,m_Fields.szFeature,
			&m_Fields.Start,&m_Fields.End,&m_Fields.ScoreOfs,szScore,&m_Fields.Strand,&cFrame,&AttribListStart);

	if(Cnt != 8)
		{
		AddErrMsg("CGTFFile::ParseLine","Errors at field %d parsing input GTF line %d",Cnt,m_CurLineNum);
		return(eBSFerrParse);
		}
	
	if(!strcmp(szScore,"."))	// if no score then set appropriate flag and default score to 0
		{	
		m_Fields.bDfltScore = true;
		m_Fields.Score = 0.0f;
		}
	else
		{	
		m_Fields.bDfltScore = false;
		m_Fields.Score = (float)atof(szScore);
		}

	if(m_Fields.Strand == '.')	// if no strand then set appropriate flag and default strand to '?'
		{	
		m_Fields.bDfltStrand = true;
		m_Fields.Strand = '?';
		}
	else
		m_Fields.bDfltStrand = false;

	if(cFrame == '.')	// if no frame then set appropriate flag and default frame to 0
		{	
		m_Fields.bDfltFrame = true;
		m_Fields.Frame = 0;
		}
	else
		{	
		m_Fields.bDfltFrame = false;
		m_Fields.Frame = cFrame - '0';
		}
    if(m_Fields.szRawLine[AttribListStart] == '\t')
		m_Fields.AttribOfs = AttribListStart+1;
	return(1);
	}
m_bEOF = true;
return(0);
}

char *
CGTFFile::GetRecord(void)				// returns complete current record
{
if(m_pGTFStream == NULL || m_bEOF)
	return(NULL);
return(m_Fields.szRawLine);
}

int 
CGTFFile::GetLineNumber(void)			// get current line number
{
if(m_pGTFStream == NULL || m_bEOF)
	return(-1);
return(m_CurLineNum);
}


int 
CGTFFile::GetCurLineLen(void)			// get current line length
{
if(m_pGTFStream == NULL || m_bEOF)
	return(-1);
return((int)strlen(m_Fields.szRawLine));
}

char *
CGTFFile::GetSeqName(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return(NULL);
return(m_Fields.szSeqName);
}

char *
CGTFFile::GetSource(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return(NULL);
return(m_Fields.szSeqSource);
}

char *
CGTFFile::GetFeature(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return(NULL);
return(m_Fields.szFeature);
}

int 
CGTFFile::GetStart(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return(-1);
return(m_Fields.Start);
}

int 
CGTFFile::GetEnd(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return(-1);
return(m_Fields.End);
}

float 
CGTFFile::GetScore(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return(-1.0f);
return(m_Fields.Score);
}

char 
CGTFFile::GetStrand(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return('?');
return(m_Fields.Strand);
}

char 
CGTFFile::GetFrame(void)
{
if(m_pGTFStream == NULL || m_bEOF)
	return('-');
return(m_Fields.Frame);
}

tsGTFFields *
CGTFFile::GetFields(void)		// returns all parsed out fields for current record
{
if(m_pGTFStream == NULL || m_bEOF)
	return(NULL);
return(&m_Fields);
}


