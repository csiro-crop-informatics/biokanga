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

#include "SQLitePSL.h"

// Following database schema is utilised
// Tables
//	TblExprs				One row for each experiment
//  TblBlatAlignment		One row for each Blat alignment
//  TblBlatAlignmentBlock	One row for each alignment block
//  TblAlignSummaries		One row for each query or target alignment summary 

// In each table the following columns are defined
//	TblExprs		One row for each experiment
//     ExprID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this experiment instance
//     ExprName VARCHAR(50) UNIQUE,			-- experiment name
//     PSLFile VARCHAR(200),				-- Blat'd query sequences were from this file
//     QueryFile VARCHAR(200),				-- Blat'd query sequences were from this file
//     TargetFile VARCHAR(200),				-- Blat'd targeted sequences were from this file
//     ExprDescr VARCHAR(1000) DEFAULT 'N/A',	-- describes experiment
//     BlatParams VARCHAR(200) DEFAULT 'N/A'   -- Blat parameters used
//     ExprType INTEGER DEFAULT 0,          -- type, currently defaulted to 0

//	TblAlignments		One row for each Blat alignment
//     AlignmentID INTEGER PRIMARY KEY ASC,	-- Uniquely identifies this alignment instance
//     ExprID INTEGER,						-- alignment was in this experiment
//     Score INTEGER,						-- Alignment score (using Blat pslScore() function)
//     Identity INTEGER,                    -- Alignment identity (using Blat 100.0 - pslCalcMilliBad(psl, TRUE) * 0.1)
//	   Matches INTEGER,						-- number of matches which aren't repeats
//     Mismatches INTEGER,					-- number of bases which do not match
//     RepMatches INTEGER,					-- number of bases which match but are also repeats
//	   NCount INTEGER,						-- umber of N bases
//	   QNumInDels INTEGER,					-- number of InDel seqs in query
//	   QBasesInDels INTEGER,				-- number of bases total in all InDels in query
//	   TNumInDels INTEGER,					-- number of InDel seqs in target
//	   TBasesInDels INTEGER,				-- number of bases total in all InDels in target
//     Strand VARCHAR(2),					-- '+' or '-' for query strand, optionally followed by target genomic strand
//     QName VARCHAR(80),					-- query sequence name
//	   QSize INTEGER,					    -- query sequence size
//	   QStart INTEGER,						-- alignment start psn in query
//     QEnd INTEGER,						-- alignment end psn in query
//     TName VARCHAR(80),					-- target sequence name
//	   TSize INTEGER,						-- target sequence size
//	   TStart INTEGER,					    -- alignment start psn in target
//	   TEnd INTEGER,					    -- alignment end psn in target
//     NumBlocks INTEGER                    -- number of alignment blocks


//	TblAlignmentBlocks		One row for each Blat alignment blocks, one or more for each alignment
//     BlockID INTEGER PRIMARY KEY ASC,		-- Uniquely identifies this block instance
//     ExprID INTEGER,						-- alignment was in this experiment
//     AlignmentID INTEGER,				-- block is in this alignment
//     BlockSize INTEGER,					-- block size
//	   QStart INTEGER,						-- starting psn of block in querry
//	   TStart INTEGER,						-- starting psn of block in target

//	TblAlignSummaries							-- One row for each target or query alignment summary
//     AlignSummaryID INTEGER PRIMARY KEY ASC,	 Uniquely identifies this alignment summary instance
//     ExprID INTEGER,						-- alignment was in this experiment
//     IsQuery INTEGER,						-- 0 if summary for target, 1 if summary for query
//     SeqName VARCHAR(80),					-- query or target sequence name
//     SeqLen INTEGER,						-- query or target sequence is this length
//	   NumAlignments INTEGER				-- if target then total number of alignments to this target, if query then total number of alignments from this query


tsStmSQL CSQLitePSL::m_StmSQL[4] = {   // 4 tables; TblExprs, TblAlignments, TblAlignmentBlocks and TblAlignSummaries
	{(char *)"TblExprs",
		(char *)"CREATE TABLE TblExprs (ExprID INTEGER PRIMARY KEY ASC, ExprName VARCHAR(50) UNIQUE,PSLFile VARCHAR(200),QueryFile VARCHAR(200), TargetFile VARCHAR(200), ExprDescr VARCHAR(1000) DEFAULT 'N/A',BlatParams VARCHAR(1000) DEFAULT 'N/A',ExprType INTEGER DEFAULT 0)",
		(char *)"INSERT INTO TblExprs (ExprName,PSLFile,QueryFile,TargetFile,ExprDescr,BlatParams,ExprType) VALUES(?,?,?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblExprs_ExprName';CREATE INDEX IF NOT EXISTS 'TblExprs_ExprName' ON 'TblExprs' ('ExprName' ASC)",
		NULL },
	{(char *)"TblAlignments",
		(char *)"CREATE TABLE TblAlignments (AlignmentID INTEGER PRIMARY KEY ASC,ExprID INTEGER,Score INTEGER,Identity INTEGER,Matches INTEGER,Mismatches INTEGER,	RepMatches INTEGER, NCount INTEGER, QNumInDels INTEGER, QBasesInDels INTEGER, TNumInDels INTEGER, TBasesInDels INTEGER,	Strand VARCHAR(2),  QName VARCHAR(80), QSize INTEGER, QStart INTEGER, QEnd INTEGER, TName VARCHAR(80), TSize INTEGER, TStart INTEGER, TEnd INTEGER,NumBlocks INTEGER, FOREIGN KEY (ExprID) REFERENCES TblExprs(ExprID))",
		(char *)"INSERT INTO TblAlignments (ExprID,Score,Identity,Matches,Mismatches,RepMatches,NCount,QNumInDels,QBasesInDels,TNumInDels,TBasesInDels,Strand,QName,QSize,QStart,QEnd,TName,TSize,TStart,TEnd,NumBlocks) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblAlignments_ExprID' ON 'TblAlignments' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_QName' ON 'TblAlignments' ('QName' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_TName' ON 'TblAlignments' ('TName' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Score' ON 'TblAlignments' ('Score' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Identity' ON 'TblAlignments' ('Identity' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Matches' ON 'TblAlignments' ('Matches' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblAlignments_ExprID' ON 'TblAlignments' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_QName' ON 'TblAlignments' ('QName' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_TName' ON 'TblAlignments' ('TName' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Score' ON 'TblAlignments' ('Score' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Identity' ON 'TblAlignments' ('Identity' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Matches' ON 'TblAlignments' ('Matches' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblAlignments_ExprID';DROP INDEX IF EXISTS 'TblAlignments_QName';DROP INDEX IF EXISTS 'TblAlignments_TName';DROP INDEX IF EXISTS 'TblAlignments_Score';DROP INDEX IF EXISTS 'TblAlignments_Identity';DROP INDEX IF EXISTS 'TblAlignments_Matches';CREATE INDEX IF NOT EXISTS 'TblAlignments_ExprID' ON 'TblAlignments' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_QName' ON 'TblAlignments' ('QName' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_TName' ON 'TblAlignments' ('TName' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Score' ON 'TblAlignments' ('Score' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Identity' ON 'TblAlignments' ('Identity' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignments_Matches' ON 'TblAlignments' ('Matches' ASC)",
		NULL },
	{(char *)"TblAlignmentBlocks",
		(char *)"CREATE TABLE TblAlignmentBlocks (BlockID INTEGER PRIMARY KEY ASC,ExprID INTEGER,AlignmentID INTEGER,BlockSize INTEGER,QStart INTEGER,TStart INTEGER, FOREIGN KEY (AlignmentID) REFERENCES TblAlignments(AlignmentID))",
		(char *)"INSERT INTO TblAlignmentBlocks (ExprID,AlignmentID,BlockSize,QStart,TStart) VALUES(?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblAlignmentBlocks_ExprID' ON 'TblAlignmentBlocks' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignmentBlocks_AlignmentID' ON 'TblAlignmentBlocks' ('AlignmentID' ASC)",
		NULL,
		(char *)"DROP INDEX IF EXISTS 'TblAlignmentBlocks_ExprID';DROP INDEX IF EXISTS 'TblAlignmentBlocks_AlignmentID';CREATE INDEX IF NOT EXISTS 'TblAlignmentBlocks_ExprID' ON 'TblAlignmentBlocks' ('ExprID' ASC); CREATE INDEX IF NOT EXISTS 'TblAlignmentBlocks_AlignmentID' ON 'TblAlignmentBlocks' ('AlignmentID' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblAlignmentBlocks_ExprID';DROP INDEX IF EXISTS 'TblAlignmentBlocks_AlignmentID'"},

	{(char *)"TblAlignSummaries",
		(char *)"CREATE TABLE TblAlignSummaries (AlignSummaryID INTEGER PRIMARY KEY ASC, ExprID INTEGER, IsQuery INTEGER, SeqName VARCHAR(80),SeqLen INTEGER,NumAlignments INTEGER)",
		(char *)"INSERT INTO TblAlignSummaries (ExprID,IsQuery,SeqName,SeqLen,NumAlignments) VALUES(?,?,?,?,?)",
		NULL,
		(char *)"CREATE INDEX IF NOT EXISTS 'TblAlignSummaries_ExprID' ON 'TblAlignSummaries' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignSummaries_SeqName' ON 'TblAlignSummaries' ('SeqName' ASC)",
		(char *)"CREATE INDEX IF NOT EXISTS 'TblAlignSummaries_ExprID' ON 'TblAlignSummaries' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignSummaries_SeqName' ON 'TblAlignSummaries' ('SeqName' ASC)",
		(char *)"DROP INDEX IF EXISTS 'TblAlignSummaries_ExprID';DROP INDEX IF EXISTS 'TblAlignSummaries_SeqName'; CREATE INDEX IF NOT EXISTS 'TblAlignSummaries_ExprID' ON 'TblAlignSummaries' ('ExprID' ASC);CREATE INDEX IF NOT EXISTS 'TblAlignSummaries_SeqName' ON 'TblAlignSummaries' ('SeqName' ASC)",
		NULL
		 }
	};


CSQLitePSL::CSQLitePSL(void)
{
m_pDB = NULL;
m_pInBuffer = NULL;		
m_pszPSLLineBuff = NULL;
m_hPSLinFile = -1;

m_NumAlignSummaries = 0;
m_UsedAlignmentSummariesSize = 0;
m_allocAlignmentSummariesSize = 0;
m_pAlignmentSummaries = NULL;
memset(m_AlignSummaryInstancesOfs,0,sizeof(m_AlignSummaryInstancesOfs));
m_NumAlignments = 0;
m_NumBlatHitsParsed = 0;
m_NumBlatHitsAccepted = 0;
m_NumBlocks = 0;
m_NumInBuffer = 0;	
m_InBuffIdx = 0;	
m_PushedBack = 0;	
m_CurLineLen = 0;
m_AllocdInBuffer = 0;
m_AllocdPSLLineBuffer = 0;
}


CSQLitePSL::~CSQLitePSL(void)
{
if(m_pDB != NULL)
	{
	sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	}
if(m_pInBuffer != NULL)
	delete m_pInBuffer;
if(m_pszPSLLineBuff != NULL)
	delete m_pszPSLLineBuff;
if(m_pAlignmentSummaries != NULL)
	{
#ifdef _WIN32
	free((UINT8 *)m_pAlignmentSummaries);
#else
	if(m_pAlignmentSummaries != MAP_FAILED)
		munmap((UINT8 *)m_pAlignmentSummaries,m_allocAlignmentSummariesSize);
#endif
	}
}

char *
CSQLitePSL::RemoveQuotes(char *pszRawText)
{
char *pSrcChr;
char *pDstChr;
bool bInSpace;
char Chr;
CUtility::TrimQuotedWhitespcExtd(pszRawText);
pSrcChr = pszRawText;
pDstChr = pSrcChr;
bInSpace = false;
while((Chr = *pSrcChr++)!= '\0')
	{
	if(Chr == '\'' || Chr == '"')
		continue;
	if(Chr == ' ' || Chr == '\t')
		{
		if(bInSpace)
			continue;
		bInSpace = true;
		}
	else
		bInSpace = false;
	*pDstChr++ = Chr;
	}
*pDstChr = '\0';
return(pszRawText);
}

sqlite3 *
CSQLitePSL::CreateDatabase(char *pszDatabase,		// database to open/create 
						bool bAppend)				// true to append onto any existing database
{
tsStmSQL *pStms;
int TblIdx;
int sqlite_error;
// note if database already exists in case bReplace is requested
struct stat TargStat;
int StatRslt = stat(pszDatabase,&TargStat);
if(StatRslt >= 0 && !bAppend)
	{
	StatRslt = -1;
	remove(pszDatabase);
	}

// try opening/creating the database
if((sqlite_error = sqlite3_open_v2(pszDatabase, &m_pDB,SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't open database: %s", sqlite3_errmsg(m_pDB));
	sqlite3_shutdown();
	m_pDB = NULL;
	return(NULL);
	}


// if required then create all tables
if(StatRslt < 0)
	{
	pStms = m_StmSQL;
	for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
		{
		pStms->pPrepInsert = NULL;
		if(pStms->pszCreateTbl == NULL)
			continue;
		if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszCreateTbl,0,0,0))!=SQLITE_OK)
			{
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
			gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszCreateTbl);   
			sqlite3_close_v2(m_pDB);
			sqlite3_shutdown();
			m_pDB = NULL;
			return(NULL);
			}
		}
	}


pStms = m_StmSQL;
for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
	{
	if(pStms->pszOpenCreateIndexes == NULL)
		continue;
	if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszOpenCreateIndexes,0,0,0))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create indexes on table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszOpenCreateIndexes);   
		sqlite3_close_v2(m_pDB);
		sqlite3_shutdown();
		m_pDB = NULL;
		return(NULL);
		}
	}

// prepare all insert statements
pStms = m_StmSQL;
for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
	{
	if(pStms->pszInsert == NULL)
		{
		pStms->pPrepInsert = NULL;
		continue;
		}
	if((sqlite_error = sqlite3_prepare_v2(m_pDB,pStms->pszInsert,-1,&pStms->pPrepInsert,NULL))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't prepare insert statement on table %s: %s", pStms->pTblName, sqlite3_errmsg(m_pDB));
		while(TblIdx > 0)
			{
			pStms -= 1;
			if(pStms->pPrepInsert != NULL)
				{
				sqlite3_finalize(pStms->pPrepInsert);
				pStms->pPrepInsert = NULL;
				}
			}
		sqlite3_close_v2(m_pDB);
		sqlite3_shutdown();
		m_pDB = NULL;
		return(NULL);
		}
	}

// initialisation for alignment sequences summaries
if(m_pAlignmentSummaries != NULL)
	{
#ifdef _WIN32
	free((UINT8 *)m_pAlignmentSummaries);
#else
	if(m_pAlignmentSummaries != MAP_FAILED)
		munmap((UINT8 *)m_pAlignmentSummaries,m_allocAlignmentSummariesSize);
#endif
	}
m_pAlignmentSummaries = NULL;
m_allocAlignmentSummariesSize = 0;
m_NumAlignSummaries = 0;
m_UsedAlignmentSummariesSize = 0;
memset(m_AlignSummaryInstancesOfs,0,sizeof(m_AlignSummaryInstancesOfs));

size_t memreq;
memreq = ((size_t)cAllocAlignSummaryInsts * (sizeof(tsAlignSummary) + (size_t)cMaxSeqNameLen));
#ifdef _WIN32
m_pAlignmentSummaries = (tsAlignSummary *) malloc(memreq);	// initial and perhaps the only allocation
if(m_pAlignmentSummaries == NULL)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateDatabase: Memory allocation of %lld bytes for alignment summary instances failed - %s",(INT64)memreq,strerror(errno));
	sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	return(NULL);
	}
#else
	// gnu malloc is still in the 32bit world and can't handle more than 2GB allocations
	m_pAlignmentSummaries = (tsAlignSummary *)mmap(NULL,memreq, PROT_READ |  PROT_WRITE,MAP_PRIVATE | MAP_ANONYMOUS, -1,0);
	if(m_pAlignmentSummaries == MAP_FAILED)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"CreateDatabase: Memory allocation of %lld bytes for alignment summary instances failed - %s",(INT64)memreq,strerror(errno));
		m_pAlignmentSummaries = NULL;
		sqlite3_close_v2(m_pDB);
		sqlite3_shutdown();
		m_pDB = NULL;
		return(NULL);
		}
#endif
m_allocAlignmentSummariesSize = memreq;
memset(m_pAlignmentSummaries,0,memreq);
return(m_pDB);
}

int
CSQLitePSL::CloseDatabase(bool bNoIndexes)
{
int TblIdx;
int Rslt = 0;
tsStmSQL *pStms;
pStms = m_StmSQL;
if(m_pDB != NULL)
	{
	if(!bNoIndexes)
		{
		for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
			{
			if(pStms->pPrepInsert == NULL)
				continue;
			sqlite3_finalize(pStms->pPrepInsert);
			pStms->pPrepInsert = NULL;
			}
		}
	Rslt = sqlite3_close_v2(m_pDB);
	sqlite3_shutdown();
	m_pDB = NULL;
	}

if(m_pAlignmentSummaries != NULL)
	{
#ifdef _WIN32
	free((UINT8 *)m_pAlignmentSummaries);
#else
	if(m_pAlignmentSummaries != MAP_FAILED)
		munmap((UINT8 *)m_pAlignmentSummaries,m_allocAlignmentSummariesSize);
#endif
	}
m_pAlignmentSummaries = NULL;
m_allocAlignmentSummariesSize = 0;
m_NumAlignSummaries = 0;
m_UsedAlignmentSummariesSize = 0;
memset(m_AlignSummaryInstancesOfs,0,sizeof(m_AlignSummaryInstancesOfs));
return(Rslt);
}

// callbacks from sqlite3_exec used to return identifier
int CSQLitePSL::ExecCallbackID(void *pCallP1, // callback function processing identifier (4th arg to sqlite3_exec())
					int NumCols,			// number of result columns 
					char **ppColValues,		// array of ptrs to column values 
					char **ppColName)		// array of ptrs to column names
{
int ValChars;
int *pID;
char *ppEnd;

// some basic validation of call back parameter values
if(pCallP1 == NULL || NumCols != 1 || ppColValues == NULL || ppColValues[0] == NULL || *ppColValues[0] == '\0')
	return(1);

pID = (int *)pCallP1;
ValChars = (int)strlen(ppColValues[0]);
*pID = strtol(ppColValues[0],&ppEnd,10);
return(0);
}


int												// errors if < eBSFSuccess, if positive then the ExprID
CSQLitePSL::CreateExperiment(char *pszExprName,	// experiment name
				char *pszPSLFile,				// alignments were parsed from this BLAT generated PSL file
				char *pszQueryFile,				// Blat'd query sequences in this file
				char *pszTargetFile,			// against targeted sequences in this file
				char *pszExprDescr,				// describes experiment
				char *pszBlatParams,			// Blat parameters used
				int ExprType)					// experiment type, currently just a place holder and defaults to 0
{
int sqlite_error;
int ExprID;

const char *pszNA = "N/A";
tsStmSQL *pStm;

if(m_pDB == NULL)
	return(eBSFerrInternal);

if(pszExprName == NULL || pszExprName[0] == '\0' ||
	pszPSLFile == NULL || pszPSLFile[0] == '\0' ||
	pszQueryFile == NULL || pszQueryFile[0] == '\0' ||
	pszTargetFile == NULL || pszTargetFile[0] == '\0')
	return(eBSFerrParams);

if(pszExprDescr == NULL || pszExprDescr[0] == '\0')
	pszExprDescr = (char *)pszNA;
if(pszBlatParams == NULL || pszBlatParams[0] == '\0')
	pszBlatParams = (char *)pszNA;

pStm = &m_StmSQL[0];

if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 1, pszExprName,(int)strlen(pszExprName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 2, pszPSLFile,(int)strlen(pszPSLFile)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, pszQueryFile,(int)strlen(pszQueryFile)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 4, pszTargetFile,(int)strlen(pszTargetFile)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 5, pszExprDescr,(int)strlen(pszExprDescr)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 6, pszBlatParams,(int)strlen(pszBlatParams)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 7, ExprType))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

	// find out the ExprID assigned to this experiment
ExprID = (int)sqlite3_last_insert_rowid(m_pDB);

return(ExprID);
}

INT32			// 20bit instance hash over the combination of parameterisation values passed into this function; if < 0 then hashing error  
CSQLitePSL::GenSummaryInstanceHash(INT32 ExprID,		// alignment summary is for alignment in this experiment
						bool bIsQuery,	// false if SeqName is for a target sequence, true if for a query sequence
						char *pszSeqName,	// NULL terminated sequence name
						UINT32  SeqLen)	// query or target sequence length
{
INT32 Hash;
char SeqNameChr;
if(ExprID < 0 || pszSeqName == NULL || pszSeqName[0] == '\0' || SeqLen == 0)
	return(-1);

Hash = ((INT32)SeqLen + ExprID) ^ (bIsQuery ? 0x05 : 0x09);
Hash &= 0x0fffff;
while((SeqNameChr = tolower(*pszSeqName++)) != '\0')
	{
	Hash <<= 3;
	Hash += SeqNameChr & 0x01f;
	Hash ^= (Hash >> 20) & 0x05;
	}
return(Hash & 0x0fffff);
}

// LocateAlignSummary locates, or creates if not already existing, a summary instance for specified parameters
// Note that the returned alignment instance will be unique for the parameterisation combination
tsAlignSummary *
CSQLitePSL::LocateAlignSummary(INT32 ExprID,	// alignment summary is for alignment in this experiment
						bool bIsQuery,			// false if SeqName is for a target sequence, true if for a query sequence
						char *pszSeqName,		// locate pre-existing, or allocate new if not pre-existing, alignment summary for bIsQuery type
						UINT32  SeqLen)			// query or target sequence length
{
int SeqNameLen;
tsAlignSummary *pAlignSummary;
INT64 AlignSummaryInstancesOfs;
INT64 NxtHashedSummaryInstOfs;
INT32 Hash;
Hash = GenSummaryInstanceHash(ExprID,bIsQuery,pszSeqName,SeqLen);
if(Hash < 0 || Hash > 0x0fffff)
	return(NULL);
SeqNameLen = (int)strlen(pszSeqName);
AlignSummaryInstancesOfs = m_AlignSummaryInstancesOfs[Hash];
if(AlignSummaryInstancesOfs > 0)		// will be 0 if no summary instance previously added with same hash
	{
	NxtHashedSummaryInstOfs = AlignSummaryInstancesOfs;
	do {
		pAlignSummary = (tsAlignSummary *)((UINT8 *)m_pAlignmentSummaries + NxtHashedSummaryInstOfs-1);
		if(pAlignSummary->ExprID == ExprID && (bool)pAlignSummary->FlgIsQuery == bIsQuery && pAlignSummary->SeqLen == SeqLen && pAlignSummary->SeqNameLen == SeqNameLen &&  !strcmp(pszSeqName,(char *)pAlignSummary->SeqName))
			return(pAlignSummary);
		}
	while((NxtHashedSummaryInstOfs = pAlignSummary->NxtHashedSummaryInstOfs) > 0);

	}

// couldn't locate existing summary instance, allocate and initialise new instance
pAlignSummary = (tsAlignSummary *)((UINT8 *)m_pAlignmentSummaries + m_UsedAlignmentSummariesSize);
memset(pAlignSummary,0,sizeof(tsAlignSummary));
pAlignSummary->ExprID = ExprID;
pAlignSummary->FlgIsQuery = bIsQuery ? 1 : 0;
pAlignSummary->HashSummaryInst = Hash;
pAlignSummary->SeqNameLen = SeqNameLen;
pAlignSummary->NumAlignments = 0;
pAlignSummary->SeqLen = SeqLen;
strcpy((char *)pAlignSummary->SeqName,pszSeqName);
pAlignSummary->AlignSummarySize = (UINT16)(sizeof(tsAlignSummary) + SeqNameLen);
pAlignSummary->NxtHashedSummaryInstOfs = AlignSummaryInstancesOfs;
m_AlignSummaryInstancesOfs[Hash] = m_UsedAlignmentSummariesSize + 1;
m_UsedAlignmentSummariesSize += pAlignSummary->AlignSummarySize;
m_NumAlignSummaries += 1;
return(pAlignSummary);
}

// AddAlignSummary locates an existing alignment summary for both the query and target (creates and initialises if not already known)
// If adding/creating for both query and target then increments NumAlignments for both
int
CSQLitePSL::AddAlignSummary(int ExprID,	// alignment summary is for this experiment
				char *pszQName,			// query sequence name, NULL or '\0' if unknown
				UINT32  QSize,			// query sequence size, 0 if unknown
				char *pszTName,			// target sequence name, NULL or '\0' if unknown
				UINT32  TSize)			// target sequence size
{
tsAlignSummary *pAQuery;
tsAlignSummary *pATarg;

// always ensure that sufficent memory has been allocated for at least 2 more summary instances to be allocated
// this is to ensure both pAQuery and pATarg will be consistent even if one or both required a alignment summary instance to be created
if((m_allocAlignmentSummariesSize - m_UsedAlignmentSummariesSize) < (2 * (sizeof(tsAlignSummary) + (size_t)(cMaxSeqNameLen))))
	{
	size_t memreq = m_allocAlignmentSummariesSize + (cAllocAlignSummaryInsts * (sizeof(tsAlignSummary) + (size_t)(cMaxSeqNameLen)));
 #ifdef _WIN32
	m_pAlignmentSummaries = (tsAlignSummary *) realloc((UINT8 *)m_pAlignmentSummaries,memreq);
#else
	m_pAlignmentSummaries = (tsAlignSummary *)mremap((UINT8 *)m_pAlignmentSummaries,m_allocAlignmentSummariesSize,memreq,MREMAP_MAYMOVE);
	if(m_pAlignmentSummaries == MAP_FAILED)
		m_pAlignmentSummaries = NULL;
#endif
	if(m_pAlignmentSummaries == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"AddAlignSummary: Memory reallocation to %lld bytes for alignment summary instances failed - %s",(INT64)memreq,strerror(errno));
		return(eBSFerrMem);
		}
	m_allocAlignmentSummariesSize = memreq;
	}

if(pszQName != NULL && pszQName[0] != '\0' && QSize != 0)
	pAQuery = LocateAlignSummary(ExprID,true,pszQName,QSize);
else
	pAQuery = NULL;

if(pszTName != NULL && pszTName[0] != '\0' && TSize != 0)
	pATarg = LocateAlignSummary(ExprID,false,pszTName,TSize);
else
	pATarg = NULL;

if(pAQuery != NULL && pATarg != NULL)
	{
	pAQuery->NumAlignments += 1;
	pATarg->NumAlignments += 1;
	}
return(eBSFSuccess);
}


int
CSQLitePSL::AddSummaryInstances2SQLite(void)
{
UINT32 Idx;
int sqlite_error;
tsStmSQL *pStm;
tsAlignSummary *pCurSummary;

if(m_pDB == NULL)
	return(eBSFerrInternal);

pCurSummary = m_pAlignmentSummaries;

for(Idx = 0; Idx < m_NumAlignSummaries; Idx++)
	{
	// add this alignment block instances
	pStm = &m_StmSQL[3];
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, pCurSummary->ExprID))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2,  pCurSummary->FlgIsQuery ? 1 : 0))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 3, (char *)pCurSummary->SeqName,(int)strlen((char *)pCurSummary->SeqName)+1,SQLITE_STATIC))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 4, pCurSummary->SeqLen))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 5, pCurSummary->NumAlignments))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	sqlite3_reset(pStm->pPrepInsert);
	pCurSummary = (tsAlignSummary *)((UINT8 *)pCurSummary + pCurSummary->AlignSummarySize);
	}
return(m_NumAlignSummaries);
}

int										// returned sequence identifier for sequence
CSQLitePSL::AddAlignment(int ExprID,		// alignment was in this experiment
					int Score,				// Alignment score (using Blat pslScore() function)
					int Identity,           // Alignment identity (using Blat 100.0 - pslCalcMilliBad(psl, TRUE) * 0.1)
					int Matches,			// number of matches which aren't repeats
					int Mismatches,			// number of bases which do not match
					int RepMatches,			// number of bases which match but are also repeats
					int NCount,				// number of N bases
					int QNumInDels,			// number of InDel seqs in query
					int QBasesInDels,		// number of bases total in all InDels in query
					int TNumInDels,			// number of InDel seqs in target
					int TBasesInDels,		// number of bases total in all InDels in target
					char *pszStrand,		// '+' or '-' for query strand, optionally followed by '+' or '-' for target genomic strand when using translated alignments
					char *pszQName,			// query sequence name
					int  QSize,				// query sequence size
					int  QStart,			// alignment start psn in query
					int  QEnd,				// alignment end psn in query
					char *pszTName,			// target sequence name
					int  TSize,				// target sequence size
					int  TStart,			// alignment start psn in target
					int  TEnd,				// alignment end psn in target
					int  NumBlocks,		// number of blocks in the alignment
					int  *pBlockSizes,		// array of sizes of each block
					int  *pQStarts,			// starting psn of each block in query
					int  *pTStarts)			// starting psn of each block in target
{
int sqlite_error;
tsStmSQL *pStm;
int AlignmentID;
int Idx;
int ChkExprID;
char szSeqTarg[200];

if(m_pDB == NULL)
	return(eBSFerrInternal);

ChkExprID = -1;

	// not a recently accessed experiment so need to check if already known to SQLite
sprintf(szSeqTarg,"select ExprID from TblExprs where ExprID = %d",ExprID);
sqlite3_exec(m_pDB,szSeqTarg,ExecCallbackID,&ChkExprID,NULL);

if(ChkExprID == -1)	// will be -1 if not already in database, treat as error
	{
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

// validated that the experiment instance identifier is already known to SQLite so can add this alignment instance
AddAlignSummary(ExprID,pszQName,QSize,pszTName,TSize);

pStm = &m_StmSQL[1];
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, Score))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, Identity))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 4, Matches))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 5, Mismatches))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 6, RepMatches))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 7, NCount))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 8, QNumInDels))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 9, QBasesInDels))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 10, TNumInDels))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 11, TBasesInDels))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 12, pszStrand,(int)strlen(pszStrand)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 13, pszQName,(int)strlen(pszQName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 14, QSize))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 15, QStart))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 16, QEnd))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_text(pStm->pPrepInsert, 17, pszTName,(int)strlen(pszTName)+1,SQLITE_STATIC))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 18, TSize))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 19, TStart))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 20, TEnd))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 21, NumBlocks))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
sqlite3_reset(pStm->pPrepInsert);

	// find out the AlignmentID assigned to this alignment instance
AlignmentID = (int)sqlite3_last_insert_rowid(m_pDB);


for(Idx = 0; Idx < NumBlocks; Idx++, pBlockSizes++,pQStarts++,pTStarts++)
	{
	// add this alignment block instances
	pStm = &m_StmSQL[2];
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 1, ExprID))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 2, AlignmentID))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 3, *pBlockSizes))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 4, *pQStarts))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_bind_int(pStm->pPrepInsert, 5, *pTStarts))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - bind prepared statement: %s", sqlite3_errmsg(m_pDB)); 
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	if((sqlite_error = sqlite3_step(pStm->pPrepInsert))!=SQLITE_DONE)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - step prepared statement: %s", sqlite3_errmsg(m_pDB));   
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	sqlite3_reset(pStm->pPrepInsert);
	}

return(AlignmentID);
}


int
CSQLitePSL::BeginPopulatingTables(void)
{
int sqlite_error;
char *pszBeginTransaction = (char *)"BEGIN TRANSACTION";
char *pszPragmaSyncOff = (char *)"PRAGMA synchronous = OFF";

char *pszPragmaJournMem = (char *)"PRAGMA journal_mode = MEMORY";

gDiagnostics.DiagOut(eDLInfo,gszProcName,"sqlite - populating tables");

// synchronous writes off
if((sqlite_error = sqlite3_exec(m_pDB,pszPragmaSyncOff,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't turn synchronous writes off: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

// bracket inserts as a single transaction
if((sqlite_error = sqlite3_exec(m_pDB,pszBeginTransaction,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't begin transactions: %s",sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

return(eBSFSuccess);
}

int 
CSQLitePSL::EndPopulatingTables(void)
{
int sqlite_error;
char *pszEndTransaction = (char *)"END TRANSACTION";
char *pszPragmaSyncOn = (char *)"PRAGMA synchronous = ON";
AddSummaryInstances2SQLite();

	// end transaction
if((sqlite_error = sqlite3_exec(m_pDB,pszEndTransaction,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't end transactions on '%s': %s", "Markers",sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Completed populating the sqlite database");

gDiagnostics.DiagOut(eDLInfo,gszProcName,"Generating indexes ...");

tsStmSQL *pStms;
pStms = m_StmSQL;
int TblIdx;
for(TblIdx = 0; TblIdx < 4; TblIdx++,pStms++)
	{
	if(pStms->pszCreateIndexes == NULL)
		continue;
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"Creating indexes on table %s ...", pStms->pTblName);
	if((sqlite_error = sqlite3_exec(m_pDB,pStms->pszCreateIndexes,0,0,0))!=SQLITE_OK)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't create indexes on table %s : %s", pStms->pTblName,sqlite3_errmsg(m_pDB));
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - statement: %s",pStms->pszCreateIndexes);   
		CloseDatabase(true);
		return(eBSFerrInternal);
		}
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"Indexes generated");
// synchronous writes off
if((sqlite_error = sqlite3_exec(m_pDB,pszPragmaSyncOn,NULL,NULL,NULL))!=SQLITE_OK)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"sqlite - can't turn synchronous writes on: %s", sqlite3_errmsg(m_pDB)); 
	CloseDatabase(true);
	return(eBSFerrInternal);
	}

CloseDatabase();
sqlite3_shutdown();
return(eBSFSuccess);
}

int
CSQLitePSL::ProcessPSL2SQLite(int PMode,		// processing mode, 0 to delete any existing then create new SQLite, 1 to append to existing SQLite
					int MinIdentity,			// minimum required identity
					int MinScore,				// minimum required score
					int MinMatches,				// minimum required base matches
					char *pszDatabase,			// SQLite database file
					char *pszExprName,			// experiment name
					char *pszPSLFile,			// alignments were parsed from this BLAT generated PSL file
					char *pszQueryFile,			// Blat'd query sequences in this file
					char *pszTargetFile,		// against targeted sequences in this file
					char *pszExprDescr,			// describes experiment
					char *pszBlatParams,		// Blat parameters used
					int ExprType)				// experiment type, currently just a place holder and defaults to 0

{
int Rslt;
int ExprID;

sqlite3_stmt *prepstatement = NULL;

m_PMode = PMode;
m_MinIdentity = MinIdentity;
m_MinScore = MinScore;
m_MinMatches = MinMatches;

sqlite3_initialize();

if((CreateDatabase(pszDatabase))==NULL)
	{
	sqlite3_shutdown();
	return(eBSFerrInternal);
	}

if((Rslt = CreateExperiment(pszExprName,pszPSLFile,pszQueryFile, pszTargetFile,pszExprDescr,pszBlatParams,ExprType)) < 1)
	{
	CloseDatabase(true);
	return(Rslt);
	}
ExprID = Rslt;

if((Rslt=BeginPopulatingTables())!=eBSFSuccess)
	return(Rslt);


if((Rslt = ProcessPSLFile(pszPSLFile,ExprID)) < eBSFSuccess)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLFile failed: %s",pszPSLFile); 
	CloseDatabase(true);
	return(Rslt);
	}

if((Rslt=EndPopulatingTables())!=eBSFSuccess)
	return(Rslt);

gDiagnostics.DiagOut(eDLInfo,gszProcName,"SQLite database ready for use");
return(eBSFSuccess);
}


//
// ProcessPSLline
// 
int
CSQLitePSL::ProcessPSLline(int ExprID)		// parse alignment PSL line which is in this experiment
{
int AlignmentID;
int Cnt;
int Psn;
int BlockIdx;
char *pChr;

int Score;
int Identity;

int Matches;			// number of matches which aren't repeats
int Mismatches;			// number of bases which do not match
int RepMatches;			// number of bases which match but are also repeats
int NCount;				// number of N bases
int QNumInDels;			// number of InDel seqs in query
int QBasesInDels;		// number of bases total in all InDels in query
int TNumInDels;			// number of InDel seqs in target
int TBasesInDels;		// number of bases total in all InDels in target
char szStrand[10];		// '+' or '-' for query strand, optionally followed by '+' or '-' for target genomic strand when using translated alignments
char szQName[cMaxSeqNameLen+1];			// query sequence name
int  QSize;				// query sequence size
int  QStart;			// alignment start psn in query
int  QEnd;				// alignment end psn in query
char szTName[cMaxSeqNameLen+1];			// target sequence name
int  TSize;				// target sequence size
int  TStart;			// alignment start psn in target
int  TEnd;				// alignment end psn in target
int  NumBlocks;			// number of blocks in the alignment
int  BlockSizes[cMaxNumPSLblocks];		// array of sizes of each block
int  QStarts[cMaxNumPSLblocks];			// starting psn of each block in query
int  TStarts[cMaxNumPSLblocks];			// starting psn of each block in target

m_pszPSLLineBuff[m_CurLineLen] = '\0';	
pChr = (char *)m_pszPSLLineBuff;
while(isspace(*pChr))
	pChr++;
if(!isdigit(*pChr))	// assume header line if not a digit	
	return(0);

Cnt = sscanf(pChr,"%d %d %d %d %d %d %d %d %n",
				&Matches,&Mismatches,&RepMatches,&NCount,&QNumInDels,
				&QBasesInDels,&TNumInDels,&TBasesInDels,&Psn);
if(Cnt != 8)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLline: Unable to parse line '%s' from input PSL file - %s", m_pszPSLLineBuff,m_szPSLinFile);
	return(eBSFerrParse);
	}

pChr += Psn;
Cnt = sscanf(pChr,"%s %s %d %d %d %s %d %d %d %d %n",
				szStrand,szQName,&QSize,&QStart,&QEnd,
				szTName,&TSize,&TStart,&TEnd,&NumBlocks,&Psn);
if(Cnt != 10)
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLline: Unable to parse line '%s' from input PSL file - %s", m_pszPSLLineBuff,m_szPSLinFile);
	return(eBSFerrParse);
	}

pChr += Psn;
for(BlockIdx = 0; BlockIdx < NumBlocks; BlockIdx++)
	{
	if(BlockIdx != (NumBlocks - 1))
		Cnt = sscanf(pChr," %d , %n",&BlockSizes[BlockIdx],&Psn);
	else
		{
		Cnt = sscanf(pChr," %d %n",&BlockSizes[BlockIdx],&Psn);
		if(Cnt == 1 && pChr[Psn] == ',')
			Psn += 1;
		}

	if(Cnt != 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLline: Unable to parse line '%s' from input PSL file - %s", m_pszPSLLineBuff,m_szPSLinFile);
		return(eBSFerrParse);
		}
	pChr += Psn;
	}
for(BlockIdx = 0; BlockIdx < NumBlocks; BlockIdx++)	
	{
	if(BlockIdx != (NumBlocks - 1))
		Cnt = sscanf(pChr," %d , %n",&QStarts[BlockIdx],&Psn);
	else
		{
		Cnt = sscanf(pChr," %d %n",&QStarts[BlockIdx],&Psn);
		if(Cnt == 1 && pChr[Psn] == ',')
			Psn += 1;
		}
	if(Cnt != 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLline: Unable to parse line '%s' from input PSL file - %s", m_pszPSLLineBuff,m_szPSLinFile);
		return(eBSFerrParse);
		}
	pChr += Psn;
	}
for(BlockIdx = 0; BlockIdx < NumBlocks; BlockIdx++)
	{
	if(BlockIdx != (NumBlocks - 1))
		Cnt = sscanf(pChr," %d , %n",&TStarts[BlockIdx],&Psn);
	else
		Cnt = sscanf(pChr," %d",&TStarts[BlockIdx]);
	if(Cnt != 1)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLline: Unable to parse line '%s' from input PSL file - %s", m_pszPSLLineBuff,m_szPSLinFile);
		return(eBSFerrParse);
		}
	pChr += Psn;
	}

Score = pslScore(Matches,Mismatches,RepMatches,QNumInDels,TNumInDels,szStrand,TSize,TStart,TEnd,NumBlocks,BlockSizes,QStarts,TStarts);			
double pslIdent = (double)pslCalcMilliBad(Matches,Mismatches,RepMatches,QNumInDels,TNumInDels,QSize,QStart,QEnd,szStrand,TSize,TStart,TEnd,NumBlocks,BlockSizes,QStarts,TStarts,true); 
Identity = (int)(100.0 - pslIdent * 0.1);
m_NumBlatHitsParsed += 1;
if(Score < m_MinScore ||
   Identity < m_MinIdentity ||
   Matches < m_MinMatches)
	return(0);	

AlignmentID = AddAlignment(ExprID,Score,Identity,Matches,Mismatches,RepMatches,NCount,QNumInDels,QBasesInDels,TNumInDels,TBasesInDels,szStrand,szQName,QSize,QStart,QEnd,szTName,TSize,TStart,TEnd,NumBlocks,BlockSizes, QStarts,TStarts);
if(AlignmentID > 0)
	m_NumBlatHitsAccepted += 1;
return(AlignmentID);
}

int		// 0: EOF -1: error >0 chr
CSQLitePSL::GetNxtPSLChr(void)	// returns buffered PSL file char
{
int Chr;
if((Chr = m_PushedBack) > 0)
	{
	m_PushedBack = 0;
	return(Chr);
	}
if(m_InBuffIdx == -1 || m_InBuffIdx >= m_NumInBuffer)
	{
	m_NumInBuffer = read(m_hPSLinFile,m_pInBuffer,cMaxInBuffSize);
	if(m_NumInBuffer <= 0)
		{
		m_InBuffIdx = -1;
		return(m_NumInBuffer);
		}
	m_InBuffIdx = 0;
	}
return(m_pInBuffer[m_InBuffIdx++]);
}


int
CSQLitePSL::ProcessPSLFile(char *pszInPSL,			// parse and load the alignments in this Blat generated PSL file into SQLite
						    int ExprID)				// the alignments are in this experiment
{
int Rslt;

if(pszInPSL == NULL || pszInPSL[0] == '\0' || ExprID <= 0)
	return(eBSFerrParams);

strncpy(m_szPSLinFile,pszInPSL,_MAX_PATH);
m_szPSLinFile[_MAX_PATH - 1] = '\0';

if(m_pInBuffer == NULL)
	{
	if((m_pInBuffer = new UINT8 [cMaxInBuffSize]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLFile: Unable to allocate %d chars buffering for input PSL file - %s", cMaxInBuffSize,m_szPSLinFile);
		return(eBSFerrMem);
		}
	m_AllocdInBuffer = cMaxInBuffSize;
	}

if(m_pszPSLLineBuff == NULL)
	{
	if((m_pszPSLLineBuff = new UINT8 [cMaxLenPSLline]) == NULL)
		{
		gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLFile: Unable to allocate %d chars buffering for input PSL lines from file - %s", cMaxLenPSLline,m_szPSLinFile);
		delete m_pInBuffer;
		m_AllocdInBuffer = 0;
		return(eBSFerrMem);
		}
	m_AllocdPSLLineBuffer = cMaxLenPSLline;
	}


#ifdef _WIN32
if((m_hPSLinFile = open(m_szPSLinFile,_O_RDWR | _O_BINARY | _O_SEQUENTIAL))==-1)
#else
if((m_hPSLinFile = open(m_szPSLinFile,O_RDWR))==-1)
#endif
	{
	gDiagnostics.DiagOut(eDLFatal,gszProcName,"ProcessPSLFile: Unable to open input file for processing - '%s' - %s", m_szPSLinFile,strerror(errno));
	delete m_pInBuffer;
	m_pInBuffer = NULL;
	m_AllocdInBuffer = 0;
	delete m_pszPSLLineBuff;
	m_pszPSLLineBuff = NULL;
	m_AllocdPSLLineBuffer = 0;
	return(eBSFerrOpnFile);
	}
gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessPSLFile: Processing input file '%s'",m_szPSLinFile); 

m_NumInBuffer = 0;
m_InBuffIdx = 0;
m_NumBlatHitsParsed = 0;
m_NumBlatHitsAccepted = 0;
Rslt = ProcessPSL(ExprID);

if(Rslt >= 0)
	gDiagnostics.DiagOut(eDLInfo,gszProcName,"ProcessPSLFile: Parsed %d alignments and accepted %d from '%s'",m_NumBlatHitsParsed,m_NumBlatHitsAccepted, m_szPSLinFile); 
m_NumBlatHitsParsed = 0;
m_NumBlatHitsAccepted = 0;
if(m_hPSLinFile != -1)
	{
	close(m_hPSLinFile);
	m_hPSLinFile = -1;
	}
if(m_pInBuffer != NULL)
	{
	delete m_pInBuffer;
	m_pInBuffer = NULL;
	m_AllocdInBuffer = 0;
	}
if(m_pszPSLLineBuff != NULL)
	{
	delete m_pszPSLLineBuff;
	m_pszPSLLineBuff = NULL;
	m_AllocdPSLLineBuffer = 0;
	}
return(Rslt);
}


int
CSQLitePSL::ProcessPSL(int ExprID)
{
int Rslt;
bool bInWhiteSpace = false;		

m_CurLineLen = 0;
Rslt = 0;
while(Rslt >= 0 && (Rslt = GetNxtPSLChr()) > 0) {
	switch((char)Rslt) {
		case '\r':			// silently slough CR - must have been generated on windows/msdos machine
			continue;

		case '\n':				// accept linefeeds - both Linux and windows are happy with lines terminated by NL
			bInWhiteSpace = false;
			if(m_CurLineLen)
				Rslt = ProcessPSLline(ExprID);
			m_CurLineLen = 0;
			continue;
		
		default:
			if(isspace(Rslt))
				{
				if(bInWhiteSpace)	// slough multiple whitespace
					continue;
				bInWhiteSpace = true;
				}
			else
				bInWhiteSpace = false;
			if(m_CurLineLen >= m_AllocdPSLLineBuffer-1)
				{
				gDiagnostics.DiagOut(eDLWarn,gszProcName,"ProcessPSL: Overlength PSL line processed from '%s'", m_szPSLinFile);
				return(eBSFerrParse);
				}
			m_pszPSLLineBuff[m_CurLineLen++] = (char)Rslt;
			continue;
			}
	}
if(!Rslt && m_CurLineLen)
	Rslt = ProcessPSLline(ExprID);
return(Rslt);
}


// following pslCalcMilliBad(), pslIsProtein() and pslScore() have been extracted from UCSC code at http://genome.ucsc.edu/FAQ/FAQblat.html#blat4
int 
CSQLitePSL::pslCalcMilliBad(int Matches,			// number of matches which aren't repeats
			int Mismatches,			// number of bases which do not match
			int RepMatches,			// number of bases which match but are also repeats
			int QNumInDels,			// number of InDel seqs in query
			int TNumInDels,			// number of InDel seqs in target
			int  QSize,				// query sequence size
			int  QStart,			// alignment start psn in query
			int  QEnd,				// alignment end psn in query
			char *pszStrand,		// '+' or '-' for query strand, optionally followed by '+' or '-' for target genomic strand when using translated alignments
			int  TSize,				// target sequence size
			int  TStart,			// alignment start psn in target
			int  TEnd,				// alignment end psn in target
			int  NumBlocks,			// number of blocks in the alignment
			int  *pBlockSizes,		// array of sizes of each block
			int  *pQStarts,			// starting psn of each block in query
			int  *pTStarts,			// starting psn of each block in target,
			bool isMrna)			// isMrna should be set to TRUE, regardless of whether the input sequence is mRNA or protein
{
int sizeMul = pslIsProtein(pszStrand,TSize,TStart,TEnd,NumBlocks,pBlockSizes,pQStarts,pTStarts) ? 3 : 1;
int qAliSize, tAliSize, aliSize;
int milliBad = 0;
int sizeDif;
int insertFactor;
int total;

qAliSize = sizeMul * (QEnd - QStart);
tAliSize = TEnd - TStart;
aliSize = min(qAliSize, tAliSize);
if (aliSize <= 0)
    return 0;
sizeDif = qAliSize - tAliSize;
if (sizeDif < 0)
    {
    if (isMrna)
        sizeDif = 0;
    else
        sizeDif = -sizeDif;
    }
insertFactor = QNumInDels;
if (!isMrna)
    insertFactor += TNumInDels;

total = (sizeMul * (Matches + RepMatches + Mismatches));
if (total != 0)
    milliBad = (int)((1000 * (Mismatches*sizeMul + insertFactor + round(3*log(1+sizeDif)))) / (double)total);
return milliBad;
}

/* is psl a protein psl (are it's blockSizes and scores in protein space) 
*/
bool 
CSQLitePSL::pslIsProtein(char *pszStrand,		// '+' or '-' for query strand, optionally followed by '+' or '-' for target genomic strand when using translated alignments
			int  TSize,				// target sequence size
			int  TStart,			// alignment start psn in target
			int  TEnd,				// alignment end psn in target
			int  NumBlocks,			// number of blocks in the alignment
			int  *pBlockSizes,		// array of sizes of each block
			int  *pQStarts,			// starting psn of each block in query
			int  *pTStarts)			// starting psn of each block in target
{
int lastBlock = NumBlocks - 1;

return  (((pszStrand[1] == '+' ) &&
     (TStart == pQStarts[lastBlock] + 3*pBlockSizes[lastBlock])) ||
    ((pszStrand[1] == '-') &&
     (TEnd == (TSize-(pTStarts[lastBlock] + 3*pBlockSizes[lastBlock])))));
}

int								// Return score for psl
CSQLitePSL::pslScore(int Matches,			// number of matches which aren't repeats
			int Mismatches,			// number of bases which do not match
			int RepMatches,			// number of bases which match but are also repeats
			int QNumInDels,			// number of InDel seqs in query
			int TNumInDels,			// number of InDel seqs in target
			char *pszStrand,		// '+' or '-' for query strand, optionally followed by '+' or '-' for target genomic strand when using translated alignments
			int  TSize,				// target sequence size
			int  TStart,			// alignment start psn in target
			int  TEnd,				// alignment end psn in target
			int  NumBlocks,			// number of blocks in the alignment
			int  *pBlockSizes,		// array of sizes of each block
			int  *pQStarts,			// starting psn of each block in query
			int  *pTStarts)			// starting psn of each block in target
{
int sizeMul = pslIsProtein(pszStrand,TSize,TStart,TEnd,NumBlocks,pBlockSizes,pQStarts,pTStarts) ? 3 : 1;

return (sizeMul * (Matches + (RepMatches>>1)) -
         sizeMul * Mismatches - QNumInDels - TNumInDels);
}



