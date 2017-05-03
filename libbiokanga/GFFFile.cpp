/*
 * CSIRO Open Source Software License Agreement (GPLv3)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENSE for the complete license information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include "../libbiokanga/commhdrs.h"
#else
#include "../libbiokanga/commhdrs.h"
#endif

CGFFFile::CGFFFile(void)
{
m_pGFFStream = NULL;
Reset();
}

CGFFFile::~CGFFFile(void)
{
Reset();
}

int 
CGFFFile::Open(char *pszFileName)		// file to open
{
Reset();

if(pszFileName == NULL || *pszFileName == '\0')
	return(eBSFerrParams);
if((m_pGFFStream = fopen(pszFileName,"r"))==NULL)
	{
	AddErrMsg("CGFFFile::Open","Unable to open GFF format file %s error: %s",pszFileName,strerror(errno));
	return(eBSFerrOpnFile);
	}
strcpy(m_szInFile,pszFileName);
memset(&m_Fields,0,sizeof(m_Fields));
m_bEOF = false;
return(eBSFSuccess);
}

int 
CGFFFile::Close(void)
{
Reset();
return(eBSFSuccess);
}

void
CGFFFile::Reset(void)
{
if(m_pGFFStream != NULL)
	{
	fclose(m_pGFFStream);
	m_pGFFStream = NULL;
	}

m_szInFile[0]='\0';

memset(&m_Fields,0,sizeof(m_Fields));
m_Fields.GFFVersion = 1;
m_CurLineNum = 0;		// current line number
m_CurNumFields = 0;		// number of fields parsed from current line
m_CurMaxFieldLen = 0;	// max len of any field in current line
m_bEOF = false;			// TRUE after last record parsed
}

char *
CGFFFile::TrimWhitespace(char *pTxt)
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
CGFFFile::GetNoteValue(char *pszAttribute)
{
int TagLen;
char *pszAttribTag;
char *pszAttribVal;
static char szValue[cGFFMaxValLen];
if(pszAttribute == NULL || m_pGFFStream == NULL || m_bEOF || m_Fields.AttribOfs == 0)
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
		TagLen = cGFFMaxValLen-1;
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
// "gene","transposon","mirna","snorna","trna","pseudogene"
// This requires a check on both the feature field and also the Note attribute as well as
// tracking on the Parent attribute
//
int					// 1 if record accepted, 0 if reached EOF, < 0 if error
CGFFFile::NextRecordOfType(etGGFFGeneType GGFFGeneType)
{
int Rslt;
static char szCurID[cMaxGeneNameLen] = "\0";
char *pszNoteVal = NULL;
char *pszGFFNoteVal;
char *pszGFFtype;
char *pszTmp;
char EndOfID;
int NoteValLen;

switch(GGFFGeneType) {
	case eGGFFany:			// any GFF record
		break;
	case eGGFFgene:			// only gene related records
		pszNoteVal = (char *)"protein_coding_gene";
		pszGFFtype = (char *)"gene";
		break;
	case eGGFFtransposon:	// only transposon related records
		pszGFFtype = (char *)"transposable_element_gene";
		pszNoteVal = NULL;
		break;
	case eGGFFmiRNA:		// only miRNA related records
		pszNoteVal = (char *)"miRNA";
		pszGFFtype = (char *)"gene";
		break;
	case eGGFFsnoRNA:		// only snoRNA related records
		pszNoteVal = (char *)"snoRNA";
		pszGFFtype = (char *)"gene";
		break;
	case eGGFFtRNA:			// only tRNA related records
		pszNoteVal = (char *)"tRNA";
		pszGFFtype = (char *)"gene";
		break;
	case eGGFFpseudogene:	// only pseudogene related records
		pszNoteVal = (char *)"pseudogene";
		pszGFFtype = (char *)"pseudogene";
		break;
	case eGGFFncRNA:	// only ncRNA/other_RNA related records
		pszNoteVal = (char *)"other_RNA";
		pszGFFtype = (char *)"gene";
		break;
	default:
		AddErrMsg("CGFFFile::NextRecordOfType","Requested GFF type '%d' not supported",GGFFGeneType);
		return(eBSFerrFeature);
	}

while((Rslt = NextLine()) > 0)
	{
	if(m_Fields.bIsMetadata)
		continue;
	if(GGFFGeneType == eGGFFany)
		return(1);
	if(!stricmp(m_Fields.szFeature,pszGFFtype) || (GGFFGeneType == eGGFFtransposon && !stricmp(m_Fields.szFeature,"transposable_element")))
		{
		if(pszNoteVal!=NULL)
			{
			if((pszGFFNoteVal = GetNoteValue((char *)"Note"))==NULL)
				continue;
			if(stricmp(pszGFFNoteVal,pszNoteVal))
				continue;
			}
		if((pszGFFNoteVal = GetNoteValue((char *)"ID"))!=NULL)
			{
			strncpy(szCurID,pszGFFNoteVal,cMaxGeneNameLen);
			szCurID[cMaxGeneNameLen-1] = '\0';
			}
		else
			szCurID[0] = '\0';
		return(1);
		}
	if(szCurID[0] == '\0')
		continue;
	if((pszGFFNoteVal = GetNoteValue((char *)"Parent"))==NULL &&
		(pszGFFNoteVal = GetNoteValue((char *)"Derives_from"))==NULL)
		continue;

	// Parent value could be multivalued so need to iterate over these values
	while(*pszGFFNoteVal != '\0')
		{
		if(*pszGFFNoteVal == ',' || *pszGFFNoteVal == ' ')
			{
			pszGFFNoteVal += 1;
			continue;
			}
		NoteValLen = 0;
		pszTmp = pszGFFNoteVal;
		while((EndOfID = *pszTmp++) && EndOfID != '.' && EndOfID != ',')
			NoteValLen += 1;
		if(!strnicmp(pszGFFNoteVal,szCurID,NoteValLen))
			return(1);
		pszGFFNoteVal += NoteValLen;
		if(EndOfID == '.')
			while((EndOfID = *pszGFFNoteVal) && EndOfID != ',')
				pszGFFNoteVal += 1;
		}
	}
return(Rslt);
}

int							// 1 if record processed, 0 if reached EOF, < 0 if error
CGFFFile::NextLine(void)			// move to next line in GFF file, skips blank lines and comment lines starting with comment char
{
return(ParseLine());
}

int 
CGFFFile::ParseLine(void)	// read next line and parse
{
int Cnt;
char *pTxt;
int AttribListStart;
char szScore[50];
char cFrame;

if(m_pGFFStream == NULL)
	{
	AddErrMsg("CGFFFile::ParseLine","GFF input stream closed");
	return(eBSFerrFileClosed);
	}
if(m_bEOF)				// last record already parsed?
	return(0);
m_Fields.AttribOfs = 0;
m_Fields.bIsMetadata = false;
m_Fields.End = 0;
m_Fields.bDfltStrand = true;
m_Fields.Frame = 0;
m_Fields.bDfltScore = true;
m_Fields.Score = 0.0f;
m_Fields.Start = 0;
m_Fields.bDfltStrand = true;
m_Fields.Strand = '?';
m_Fields.szFeature[0] = '\0';
m_Fields.szRawLine[0] = '\0';
m_Fields.szSeqName[0] = '\0';
m_Fields.szSeqSource[0] = '\0';
AttribListStart = 0;

while(fgets(m_Fields.szRawLine,sizeof(m_Fields.szRawLine)-1,m_pGFFStream)!=NULL)
	{
	m_CurLineNum += 1;
	pTxt = TrimWhitespace(m_Fields.szRawLine);
	m_CurLineLen = (int)strlen(pTxt);
	if(*pTxt=='\0' || (*pTxt=='#' && (pTxt[1] != '#' || pTxt[1] == '\0')))	// simply slough lines which were just whitespace or start with '#'
		continue;
	if(pTxt[0] == '#' && pTxt[1] == '#')
		{
		m_Fields.bIsMetadata = true;
		if(!strnicmp(&pTxt[2],"gff-version",11))
			m_Fields.GFFVersion = atoi(&pTxt[14]);
		return(1);
		}
	// parse out the mandatory fields
	Cnt = sscanf(pTxt,"%34s\t%34s\t%34s\t%d\t%d\t%25s\t%c\t%c%n",
			m_Fields.szSeqName,m_Fields.szSeqSource,m_Fields.szFeature,
			&m_Fields.Start,&m_Fields.End,szScore,(char *)&m_Fields.Strand,&cFrame,&AttribListStart);

	if(Cnt != 8)
		{
		AddErrMsg("CGFFFile::ParseLine","Errors at field %d parsing input GFF line %d",Cnt,m_CurLineNum);
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
CGFFFile::GetRecord(void)				// returns complete current record
{
if(m_pGFFStream == NULL || m_bEOF)
	return(NULL);
return(m_Fields.szRawLine);
}

int
CGFFFile::GetGFFversion(void)
{
if(m_pGFFStream == NULL || m_bEOF)
	return(-1);
return(m_Fields.GFFVersion);
}

int 
CGFFFile::GetLineNumber(void)			// get current line number
{
if(m_pGFFStream == NULL || m_bEOF)
	return(-1);
return(m_CurLineNum);
}


int 
CGFFFile::GetCurLineLen(void)			// get current line length
{
if(m_pGFFStream == NULL || m_bEOF)
	return(-1);
return((int)strlen(m_Fields.szRawLine));
}

bool 
CGFFFile::IsMetadataLine(void)			// returns true if current line contains metadata
{
if(m_pGFFStream == NULL || m_bEOF)
	return(false);
return(m_Fields.bIsMetadata);
}

char *
CGFFFile::GetMetadata(void)
{
if(m_pGFFStream == NULL || m_bEOF || !m_Fields.bIsMetadata)
	return(NULL);
return(m_Fields.szRawLine);
}

char *
CGFFFile::GetSeqName(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return(NULL);
return(m_Fields.szSeqName);
}

char *
CGFFFile::GetSource(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return(NULL);
return(m_Fields.szSeqSource);
}

char *
CGFFFile::GetFeature(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return(NULL);
return(m_Fields.szFeature);
}

int 
CGFFFile::GetStart(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return(-1);
return(m_Fields.Start);
}

int 
CGFFFile::GetEnd(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return(-1);
return(m_Fields.End);
}

float 
CGFFFile::GetScore(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return(-1.0f);
return(m_Fields.Score);
}

char 
CGFFFile::GetStrand(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return('?');
return(m_Fields.Strand);
}

char 
CGFFFile::GetFrame(void)
{
if(m_pGFFStream == NULL || m_bEOF || m_Fields.bIsMetadata)
	return('-');
return(m_Fields.Frame);
}
