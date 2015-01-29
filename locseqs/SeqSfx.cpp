#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stdafx.h"

#if _WIN32

#include "../conservlib/commhdrs.h"

#else
#include "../libbiokanga/commhdrs.h"
#endif

#include "SeqSfx.h"

const int cAllocSfxEls = 1000000;			// allocate new suffix index elements in this many increments
const int cAllocConcatSeqs = 1000000;		// allocate new concatenated probe sequences in this many base increments
const int cAllocSeqEls = 10000;				// allocate for new tsSeqHdr array elements in this many increments

CSeqSfx::CSeqSfx(void)
{
m_pConcatSeqs = NULL;
m_pSfxEls = NULL;
m_pPutCores = NULL;
m_pSeqHdrs = NULL;
Reset();
}

CSeqSfx::~CSeqSfx(void)
{
Reset();
}

void
CSeqSfx::Reset(void)
{
if(m_pConcatSeqs != NULL)		// pts to concatenated target sequences
	{
	delete m_pConcatSeqs;
	m_pConcatSeqs = NULL;
	}
if(m_pSfxEls != NULL)			// pts to generated sfx array
	{
	delete m_pSfxEls;
	m_pSfxEls = NULL;
	}
if(m_pPutCores != NULL)			// pts to putative cores
	{
	delete m_pPutCores;
	m_pPutCores = NULL;
	}
if(m_pSeqHdrs != NULL)			// pts to array of seqence headers
	{
	delete m_pSeqHdrs;
	m_pSeqHdrs = NULL;
	}

m_NumSfxEls = 0;		// number of elements in m_pSfxEls;
m_ConcatSeqLens = 0;	// total concatenated target sequences length including all EOS indicators
m_AllocSeqLen = 0;		// currently allocated 
m_NumSeqs = 0;			// number of sequences currently concatenated
m_AllocSeqEls = 0;		// currently allocated
m_AllocSfxEls = 0;		// currently allocated 
m_NumPutCores = 0;		// number of putative cores
m_AllocPutCoreEls = 0;	// number currently allocd for m_pPutCores
}

int 
CSeqSfx::Add(unsigned int TargetRef,	// target to associate with any hit to this sequence
			int SeqLen,// sequence length
			etSeqBase *pSeq)// probe sequence
{
sCSfxSeqs *pNewSfx;
etSeqBase *pBases;
tsSeqHdr *pNewSeqHdr;
int Idx;

// quick parameter check
if(SeqLen < 1 || pSeq == NULL)
	return(-1);

// allocate as may be required for sequence headers
if(m_pSeqHdrs == NULL || m_NumSeqs >= m_AllocSeqEls)
	{
	int AllocEls = m_AllocSeqEls + cAllocSeqEls;
	if((pNewSeqHdr = new tsSeqHdr [AllocEls])==NULL)
		return(-1);

	if(m_pSeqHdrs != NULL)
		{
		if(m_NumSeqs)
			memcpy(pNewSeqHdr,m_pSeqHdrs,m_NumSeqs * sizeof(tsSeqHdr));
		delete m_pSeqHdrs;
		}
	else
		m_NumSeqs = 0;
	m_AllocSeqEls = AllocEls;
	m_pSeqHdrs = pNewSeqHdr;
	}

// allocate as may be required for holding target sequences concatenated with EOS markers separating each sequence
if(m_pConcatSeqs == NULL || (m_ConcatSeqLens + SeqLen + 1) >= m_AllocSeqLen)
	{
	int AllocSeqLen = m_ConcatSeqLens + SeqLen + 1 + cAllocConcatSeqs;
	if((pBases = new etSeqBase [AllocSeqLen])==NULL)
		return(-1);

	if(m_pConcatSeqs != NULL)
		{
		if(m_ConcatSeqLens)
			memcpy(pBases,m_pConcatSeqs,m_ConcatSeqLens);
		delete m_pConcatSeqs;
		}
	else
		m_ConcatSeqLens = 0;
	m_pConcatSeqs = pBases;
	m_AllocSeqLen = AllocSeqLen;
	}

// allocate as may be required for suffix array elements
if(m_pSfxEls == NULL || (m_NumSfxEls + SeqLen + 1) >= m_AllocSfxEls)
	{
	int AllocEls = m_AllocSfxEls + SeqLen + 1 + cAllocSfxEls;
	if((pNewSfx = new sCSfxSeqs [AllocEls])==NULL)
		return(-1);

	if(m_pSfxEls != NULL)
		{
		if(m_NumSfxEls)
			memcpy(pNewSfx,m_pSfxEls,m_NumSfxEls * sizeof(sCSfxSeqs));
		delete m_pSfxEls;
		}
	else
		m_NumSfxEls = 0;
	m_pSfxEls = pNewSfx;
	m_AllocSfxEls = AllocEls;
	}

// initialise header for sequence
pNewSeqHdr = &m_pSeqHdrs[m_NumSeqs++];
pNewSeqHdr->SeqIdx = m_NumSeqs;
pNewSeqHdr->TargetRef = TargetRef;
pNewSeqHdr->SeqLen = SeqLen;

// concatenate sequence on to previous sequences with additional eBaseEOS to act as sequence separator
pBases = &m_pConcatSeqs[m_ConcatSeqLens];
memcpy(pBases,pSeq,SeqLen);
pBases[SeqLen] = eBaseEOS;
m_ConcatSeqLens += SeqLen + 1;

// initialise suffix array elements with sequence header index and ptr to bases from sequence
pNewSfx = &m_pSfxEls[m_NumSfxEls];
for(Idx = 0; Idx <= SeqLen; Idx++, pNewSfx++)	// iterate <= SeqLen as sequence will have had eBaseEOS appended as a sequence separator indicator
	{
	pNewSfx->SeqIdx = m_NumSeqs;
	pNewSfx->SeqOfs = Idx;
	}
m_NumSfxEls += SeqLen + 1;

return(m_NumSeqs);
}

// GenSfx
// generate suffixes over all concatenated target sequences
int 
CSeqSfx::GenSfx(void)		
{
sCSfxSeqs *pSfxEl;
etSeqBase *pBase;
tsSeqHdr *pSeqHdr;

int Idx;
pSfxEl = m_pSfxEls;
pBase = m_pConcatSeqs;
pSeqHdr = m_pSeqHdrs;
for(Idx = 0; Idx < m_NumSfxEls; Idx++,pSfxEl++,pBase++)
	{
	if(pSfxEl->SeqOfs == 0)
		{
		pSeqHdr->pStart = pBase;
		pSeqHdr += 1;
		}
	pSfxEl->pBase = pBase;
	}
qsort(m_pSfxEls,m_NumSfxEls,sizeof(sCSfxSeqs),SortSeqSfxs);

return(m_NumSfxEls);
}


// LocateFirstExact
// Locates first instance (lowest pSfxArray[] positional index) of exactly matching probe 
int			// index in pSfxArray of exactly matching probe or -1 if no match					
CSeqSfx::LocateFirstExact(etSeqBase *pProbe,  // pts to probe sequence
				  unsigned int ProbeLen,		// probe length to exactly match over
				  int SfxLo,					// low index in pSfxArray
				  int SfxHi)					// high index in pSfxArray
{
etSeqBase *pEl1;
etSeqBase *pEl2;

etSeqBase b1;
etSeqBase b2;
int CurMatchLen;
int Max2Compare;
int NumMatches;
int TargPsn;
int LowPsn;

if((*pProbe & ~cRptMskFlg) > eBaseN)
	return(-1);

NumMatches = 0;
do {
	TargPsn = (SfxLo + SfxHi) / 2;
	pEl2 = m_pSfxEls[TargPsn].pBase;
	Max2Compare = ProbeLen; 
	CurMatchLen = 0;
	pEl1 = pProbe;

	while(Max2Compare-- && (b1 = (*pEl1++ & ~cRptMskFlg)) == (b2 = (*pEl2++ & ~cRptMskFlg)) && b1 <= eBaseT)
		CurMatchLen++;

	if(b1 > eBaseT)	// probes only allowed to contain nucleotides a,c,g,t
		return(-1);

	// manage to get minimum core length?
	if(CurMatchLen == (int)ProbeLen)
		{
		if(TargPsn == 0 || SfxLo == TargPsn)
			return(TargPsn);
		LowPsn = LocateFirstExact(pProbe,ProbeLen,SfxLo,TargPsn - 1);
		return(LowPsn < 0 ? TargPsn : LowPsn);
		}

	if(b2 > eBaseT)
		SfxHi = TargPsn-1;
	else
		{
		if(b1 < b2)
			SfxHi = TargPsn - 1;
		else
			SfxLo = TargPsn+1;
		}
	}
while(SfxHi >= SfxLo);

return(-1);	// unable to locate any instance of pProbe
}



// LocateNearExacts
// Returns count of all nearly exact maximally extended matching sequences
// These matching sequences must contain an exactly matching core of MinCoreLen and
// left/right flanking sequences which are maximally extended with at most LeftMaxMismatches/RightMaxMismatches
// The flanking sequences are extended for at most LeftMaxExtend and RightMaxExtend bases.
int						
CSeqSfx::LocateNearExacts(sAMProbe *pProbe,		// probe
 				  unsigned int ProbeLen,	    // probe length
				  etSeqBase *pProbeSeq,			// probe sequence
				  unsigned int MinCoreLen,	    // minimal core length from which to extend left + right
				  unsigned int LeftMaxExtend,   // maximum length to extend left
				  unsigned int RightMaxExtend,  // maximum length to extend right
  				  unsigned int LeftMaxMismatches,  // total mismatches allowed in left flank
   				  unsigned int RightMaxMismatches, // total mismatches allowed in right flank
				  unsigned int MinHitLen,	       // putative hits must be of at least this length
				  
				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int (* FilterCore)(sAMProbe *pProbe,				// probe
										unsigned int ProbeOfs,		// probe offset at which probe hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int ProbeCoreOfs,	// probe offset at which exact core starts
										unsigned int ProbeCoreLen,	// length of exact core
										unsigned int NumLeftMismatches, // number of mismatches in left flank
										unsigned int NumRightMismatches, // number of mismatches in right flank
										unsigned int TargetRef,		// user assigned target to associate with any hit to this sequence
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				  
				 // following is call back on all accepted hits of probe into suffix entry
				  int (* ProcessCore)(sAMProbe *pProbe,			// probe
										unsigned int ProbeOfs,		// probe offset at which hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int TargetRef,		// user assigned target to associate with any hit to this sequence
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				etRPTMasking RPTMasking,							// how to interpret probe bases with cRptMskFlg set
				unsigned int MaskPercent)						// if matches have more than this percent repeats then slough
{
int Rslt;
int TargPsn;	// index in pSfxArray of current suffix sequence being processed
unsigned int CurProbeOfs;
unsigned int CurMatchLen;
unsigned int CurLeftTargOfs;
unsigned int CurLeftProbeOfs;

unsigned int ExactMatchCoreLen;
unsigned int LeftExactCoreStart;
unsigned int RightExactCoreEnd;

unsigned int MinSeedLen;
unsigned int DeltaProbeOfs;

int NumProcessed;	// total number of cores after filtering
tsCPutCore *pCore;
bool bTerm;
int Idx;

tsSeqHdr *pSeqHdr;

etSeqBase *pPrb;
etSeqBase *pTrg;
etSeqBase *pPrbr;
etSeqBase *pTrgr;
etSeqBase *pMaxPrbr;
etSeqBase *pMaxTrgr;

etSeqBase *pSeq2Loc;


unsigned int CurProbeLen;

etSeqBase b1;
etSeqBase b2;

int Max2Compare;

unsigned int LeftMismatches;
unsigned int RightMismatches;


if(MinHitLen > ProbeLen)
	return(-1);

if(MinCoreLen <= 12)
	{
	MinSeedLen = MinCoreLen;
	DeltaProbeOfs = 1;
	}
else		// else minhitlen must be at >= 13
	{
	if(MinCoreLen <= 49)
		{
		DeltaProbeOfs = MinCoreLen / 4;
		MinSeedLen = MinCoreLen - DeltaProbeOfs;
		if(MinSeedLen < 12)	// can only be if minhitlen < 16
			{
			MinSeedLen = 12;
			DeltaProbeOfs = MinCoreLen - MinSeedLen;
			}
		}
	else	// else minhitlen must be at >= 50 so MinSeedLen will be at least 25
		{
		DeltaProbeOfs = MinCoreLen / 2;
		MinSeedLen = MinCoreLen - DeltaProbeOfs;
		}
	}

m_NumPutCores = 0;			// clears out any previous putative cores
NumProcessed = 0;

CurProbeOfs = 0;
bTerm = false;

do {
	pSeq2Loc = &pProbeSeq[CurProbeOfs];			// initial start

	// putative cores which contain any bases other than canonical (a,c,g or t) are simply sloughed
	if(FiltOutRepeats(eRPTHignor,0,pSeq2Loc,MinSeedLen))
		continue;

	// locate first core of at least MinSeedLen which exactly matches
	if((TargPsn = LocateFirstExact(pSeq2Loc,MinSeedLen,0,m_NumSfxEls-1)) < 0)
		continue;

	do  {
		CurProbeLen = MinSeedLen;
		pTrg = m_pSfxEls[TargPsn].pBase;				// pts to target base
		CurLeftTargOfs = m_pSfxEls[TargPsn].SeqOfs;		// offset in target sequence

		pPrbr = &pProbeSeq[CurProbeOfs + MinSeedLen];
		pTrgr = pTrg + MinSeedLen;

		pSeqHdr = &m_pSeqHdrs[m_pSfxEls[TargPsn].SeqIdx-1];	// pt to header for target sequence
		pMaxTrgr = pSeqHdr->pStart + pSeqHdr->SeqLen - 1; // pts to last base in probe
		pMaxPrbr = &pProbeSeq[ProbeLen-1];			// pts to last base in Seq

		// firstly try to extend left
		CurMatchLen = MinSeedLen;
		ExactMatchCoreLen = MinSeedLen;
		CurLeftProbeOfs = CurProbeOfs;
		LeftExactCoreStart = CurProbeOfs;
		RightExactCoreEnd = LeftExactCoreStart + MinSeedLen - 1;
		LeftMismatches = 0;
		if(CurLeftProbeOfs > 0 && CurLeftTargOfs > 0)
			{
			pPrb = &pProbeSeq[CurLeftProbeOfs - 1];
			pTrg -= 1;
			do	{
				b1 = (*pPrb-- & ~cRptMskFlg);
				b2 = (*pTrg-- & ~cRptMskFlg);
				if(b1 > eBaseT || b2 > eBaseT)			// can't handle non-canonical
					break;
				if((b1 != b2) && ++LeftMismatches > LeftMaxMismatches)
					{
					LeftMismatches -= 1;
					break;
					}
				if(b1 == b2 && !LeftMismatches)
					{
					ExactMatchCoreLen++;
					LeftExactCoreStart--;
					}
				CurMatchLen++;
				CurLeftProbeOfs -= 1;
				CurLeftTargOfs -= 1;
				}
			while(CurLeftProbeOfs > 0 && CurLeftTargOfs > 0);
			}

				// after extending left then try to extend right
		RightMismatches = 0;
		while(pPrbr <= pMaxPrbr && pTrgr <= pMaxTrgr)
			{
			b1 = (*pPrbr++ & ~cRptMskFlg);
			b2 = (*pTrgr++ & ~cRptMskFlg);
			if(b1 > eBaseT || b2 > eBaseT)			// can't handle non-canonical
				break;
			if((b1 != b2) && ++RightMismatches > RightMaxMismatches)
				{
				RightMismatches -= 1;
				break;
				}
			if(b1 == b2 && !RightMismatches)
				{
				ExactMatchCoreLen += 1;
				RightExactCoreEnd += 1;
				}
			CurMatchLen++;
			}

		if(ExactMatchCoreLen >= MinCoreLen && CurMatchLen >= MinHitLen)
			{
				// apply filtering for repeats and non-canonical bases
			if(!FiltOutRepeats(RPTMasking,MaskPercent,&pProbeSeq[CurLeftProbeOfs],CurMatchLen))
				{			
					// check if caller wants to filter this putative approximate match out
				if(FilterCore != NULL)
					{
					if((Rslt = (*FilterCore)(pProbe,
									CurLeftProbeOfs,
									CurMatchLen,
									CurLeftTargOfs,
									ExactMatchCoreLen,
									LeftMismatches,
									RightMismatches,
									pSeqHdr->TargetRef,
									CurLeftTargOfs,
									pSeqHdr->pStart)) < 0)
							return(Rslt);
					if(!Rslt)
						{
						bTerm = true;
						continue;
						}
					}
				else
					Rslt = 1;
				if(Rslt > 0 && (Rslt =AddPutCore(CurLeftProbeOfs,	// probe offset at which hit starts
							CurMatchLen,							// number of bases hit
							pSeqHdr->SeqIdx,						// header index
							CurLeftTargOfs)) < 1)					// Seq loci hit
					return(Rslt);
				}
			}

		// any more suffixes to check
		if(++TargPsn >= m_NumSfxEls)
			break;
		// check if probe matches next suffix
		Max2Compare = MinSeedLen;
		CurMatchLen = 0;
		pTrg = m_pSfxEls[TargPsn].pBase;
		pPrb = pSeq2Loc;
		while(Max2Compare-- && (b1 = (*pPrb++ & ~cRptMskFlg)) == (b2 = (*pTrg++ & ~cRptMskFlg)) && b1 <= eBaseT)
			CurMatchLen++;
		}
	while(!bTerm && CurMatchLen >= MinSeedLen);
	}
while(!bTerm && (CurProbeOfs += DeltaProbeOfs) <= (ProbeLen - MinSeedLen));

NumProcessed = FiltPutCores();
if(NumProcessed <= 0)
	return(NumProcessed);
pCore = m_pPutCores;
int NumAccepted = 0;
for(Idx = 0; Idx < NumProcessed; Idx++,pCore++)
	{
	if(pCore->HitLen)
		{
		if((Rslt = (*ProcessCore)(pProbe,pCore->ProbeOfs,pCore->HitLen,pSeqHdr->TargetRef,pCore->SeqOfs,pSeqHdr->pStart))<=0)
			return(Rslt);
		NumAccepted += 1;
		}
	}
return(NumAccepted);
}


int
CSeqSfx::AddPutCore(unsigned int ProbeOfs,	// probe offset at which hit starts 
	unsigned int HitLen,			// number of bases hit
	UINT32 SeqIdx,					// identifies which Seq (1..N)
	unsigned int TargOfs)			// target loci hit
{
tsCPutCore *pCore;
int CoreIdx;
int ProbeEndOfs = ProbeOfs + HitLen - 1;
int TargEndOfs = TargOfs + HitLen - 1;

// iterate cores to determine if core is contained within existing core, or contains an existing core
// if contained within an existing then don't add, if containing then replace existing 
if(m_pPutCores != NULL && m_NumPutCores)
	{
	pCore = m_pPutCores;
	for(CoreIdx = 0; CoreIdx < m_NumPutCores; CoreIdx++,pCore++)
		{
		if(pCore->SeqIdx != SeqIdx ||		// if not same target or
			pCore->HitLen == 0)			// core marked to be sloughed
			continue;					// then check next core
		if(pCore->ProbeOfs <= (UINT32)ProbeOfs && (pCore->ProbeOfs + pCore->HitLen - 1) >= (UINT32)ProbeEndOfs)
			{
			if(pCore->SeqOfs <= TargOfs && (pCore->SeqOfs + pCore->HitLen - 1) >= (UINT32)TargEndOfs)
				return(m_NumPutCores);
			}

		if(pCore->ProbeOfs >= (UINT32)ProbeOfs && (pCore->ProbeOfs + pCore->HitLen - 1) <= (UINT32)ProbeEndOfs)
			{
			if(pCore->SeqOfs >= TargOfs && (pCore->SeqOfs + pCore->HitLen - 1) <= (UINT32)TargEndOfs)
				{
				pCore->ProbeOfs = ProbeOfs;
				pCore->HitLen = HitLen;
				pCore->SeqOfs = TargOfs;
				return(m_NumPutCores);
				}
			}
		}
	}

// core not contained within any other core so add
if(m_pPutCores == NULL || m_NumPutCores >= m_AllocPutCoreEls)
	{
	if((pCore = new tsCPutCore [cAllocPutCores + m_AllocPutCoreEls])==NULL)
		return(eBSFerrMem);

	if(m_pPutCores != NULL)
		{
		memcpy(pCore,m_pPutCores,m_NumPutCores * sizeof(tsPutCore));
		delete m_pPutCores;
		m_AllocPutCoreEls += cAllocPutCores;
		}
	else
		{
		m_AllocPutCoreEls = cAllocPutCores;
		m_NumPutCores = 0;
		}
	m_pPutCores = pCore;
	}
pCore = &m_pPutCores[m_NumPutCores++];
pCore->SeqIdx = SeqIdx;
pCore->ProbeOfs = ProbeOfs;
pCore->HitLen = HitLen;
pCore->SeqOfs = TargOfs;
return(m_NumPutCores);
}


// FiltContainedPutCores
// Marks for deletion (sets core length to 0) any cores starting from NxtProbeIdx which are 
// internally contained within cores starting at CurProbeIdx
int								// returns number of cores to be deleted
CSeqSfx::FiltContainedPutCores(int CurProbeIdx, 
								 int NxtProbeIdx)
{
int Num2Mark;

unsigned int CurStartOfs;
unsigned int CurEndOfs;
unsigned int HitLen;
unsigned int CurTargOfs;
unsigned int CurTargEndOfs;

tsCPutCore *pCurCore;
tsCPutCore *pNxtCore;
pCurCore = &m_pPutCores[CurProbeIdx];

// nothing to do if probe already marked for removal?
if(!(HitLen = pCurCore->HitLen))
	return(0);

CurStartOfs = pCurCore->ProbeOfs;
CurEndOfs = CurStartOfs + pCurCore->HitLen - 1;
CurTargOfs = pCurCore->SeqOfs;
CurTargEndOfs = CurTargOfs + pCurCore->HitLen - 1;
pNxtCore = &m_pPutCores[NxtProbeIdx];
Num2Mark = 0;
do	{
	if(pCurCore->SeqIdx < pNxtCore->SeqIdx)	// if next core for sequence after current then there can't be any overlap
			break;

	if(pNxtCore->HitLen == 0 || pCurCore->SeqIdx != pNxtCore->SeqIdx)	// if target already marked or not target sequence of interest then skip to next
		{
		pNxtCore+=1;
		continue;
		}

	if(pNxtCore->ProbeOfs > CurEndOfs)	// if next starts after current ends then subsequent probes can't be contained
		break;

	if(CurStartOfs <= pNxtCore->ProbeOfs &&
		CurEndOfs >= (pNxtCore->ProbeOfs + pNxtCore->HitLen - 1) && // if target core is contained within probe core
		CurTargOfs <= pNxtCore->SeqOfs &&							// and target loci is contained
		CurTargEndOfs >= (pNxtCore->SeqOfs + pNxtCore->HitLen - 1))
		{
		pNxtCore->HitLen = 0;	// mark target core to be removed
		Num2Mark++;
		}

	pNxtCore+=1;
	}
while(++NxtProbeIdx < m_NumPutCores);

return(Num2Mark);
}

// FiltOverlapPutCores
// Marks for deletion the shorter of any cores which overlap both on the CurProbeIdx and NxtProbeIdx..m_NumPutCores
int							// returns number of cores marked for deletion
CSeqSfx::FiltOverlapPutCores(int CurProbeIdx, int NxtProbeIdx)
{
int Num2Mark;

unsigned int CurStartOfs;
unsigned int CurEndOfs;
unsigned int HitLen;
unsigned int CurTargOfs;
tsCPutCore *pCurCore;
tsCPutCore *pNxtCore;
pCurCore = &m_pPutCores[CurProbeIdx];
if(!(HitLen = pCurCore->HitLen))
	return(0);
CurStartOfs = pCurCore->ProbeOfs;
CurEndOfs = CurStartOfs + HitLen - 1;
CurTargOfs = pCurCore->SeqOfs;
pNxtCore = &m_pPutCores[NxtProbeIdx];
Num2Mark = 0;

do	{
	if(pCurCore->SeqIdx < pNxtCore->SeqIdx)	// if next core for sequence after current then there can't be any overlap
			break;

	if(pNxtCore->HitLen == 0 ||	// on to next if already marked for deletion or different sequence
		pCurCore->SeqIdx != pNxtCore->SeqIdx)
		{
		pNxtCore+=1;
		continue;
		}

	if(pNxtCore->ProbeOfs > CurEndOfs)		// if next core start is after current core ends then can't be any overlap
		break;

	if((CurStartOfs <= pNxtCore->ProbeOfs) && // already know that pNxtCore->ProbeOfs <= CurEndOfs so here we are checking for any overlap on the probe cores
		((pNxtCore->SeqOfs + pNxtCore->HitLen - 1) >= CurTargOfs) && // assuming probe core overlaps then here we are checking for overlap on the target
		(pNxtCore->SeqOfs <= (CurTargOfs + HitLen - 1)))
		{									  // have an overlap on both the probe core and target sequence loci 
		// retain longest core, shortest to be deleted
		if(pNxtCore->HitLen > HitLen)
			{
			pCurCore->HitLen = 0;
			Num2Mark++;
			break;
			}
		else
			pNxtCore->HitLen = 0;	// mark target core to be removed
		Num2Mark++;
		}
	pNxtCore+=1;
	}
while(++NxtProbeIdx < m_NumPutCores);

return(Num2Mark);
}


// FiltPutCores
// Filter cores for contained hits from probes onto target sequence, optionally also remove overlaps
int
CSeqSfx::FiltPutCores(bool bRemoveOverlaps)
{
int MarkIdx;
int NumMarked;
bool bResort;			// will be set true if any cores marked and resort is required
tsCPutCore *pCore;

if(m_NumPutCores <= 1)
	return(m_NumPutCores);

qsort(m_pPutCores,m_NumPutCores,sizeof(tsCPutCore),SortPutCores);

// mark any cores contained within other cores
bResort = false;
do {
	NumMarked = 0;
	for(MarkIdx = 0; MarkIdx < m_NumPutCores-1; MarkIdx++)
		NumMarked += FiltContainedPutCores(MarkIdx,MarkIdx+1);
	if(NumMarked)
		bResort = true;
	}
while(NumMarked > 0);

if(bRemoveOverlaps)	// are overlaps to be removed?
	{
	if(bResort)
		{
		// sort to put marked cores last and then update m_NumPutCores with number of cores not marked
		qsort(m_pPutCores,m_NumPutCores,sizeof(tsCPutCore),SortPutCores);
		pCore = m_pPutCores;
		for(MarkIdx = 0; MarkIdx < m_NumPutCores; MarkIdx++,pCore++)
			if(pCore->HitLen == 0)
				break;
		m_NumPutCores = MarkIdx;
		}

	bResort = false;
	// mark any overlap cores
	do {
		NumMarked = 0;
		for(MarkIdx = 0; MarkIdx < m_NumPutCores-1; MarkIdx++)
			NumMarked += FiltOverlapPutCores(MarkIdx,MarkIdx+1);
		if(NumMarked)
			bResort = true;
		}
	while(NumMarked > 0);
	}

if(bResort)
	{
	// sort to put marked cores last and then update m_NumPutCores with number of cores not marked
	qsort(m_pPutCores,m_NumPutCores,sizeof(tsCPutCore),SortPutCores);
	pCore = m_pPutCores;
	for(MarkIdx = 0; MarkIdx < m_NumPutCores; MarkIdx++,pCore++)
		if(pCore->HitLen == 0)
			break;
	m_NumPutCores = MarkIdx;
	}
return(m_NumPutCores);
}


// FiltOutRepeats
// Returns true if sequence is to be filtered out either because it contains a non-canonical base or the percentage of repeat masked bases is too high
bool
CSeqSfx::FiltOutRepeats(etRPTMasking RPTMasking,	// how to interpret cRptMskFlg'd bases
						 unsigned int MaskPercent,	// filter out if percentage of repeats is above this percentage (0-100)
						 etSeqBase *pSeq,			// pts to sequence
						 int SeqLen)				// sequence length
{
etSeqBase Base;
unsigned int Rpts = 0;
unsigned int Percent;
int Len = SeqLen;
switch(RPTMasking) {
	case eRPTHignor:	// treat all bases as not a repeat but still check for non-canonical bases
		while(Len--)
			if((*pSeq++ & ~cRptMskFlg)  > eBaseT)
				return(true);
		return(false);

	case eRPTHmasked:	// treat any base cRptMskFlg as a repeat masked base
		while(Len--)
			{
			if(((Base = *pSeq++) & ~cRptMskFlg)  > eBaseT)
				return(true);
			if(Base & cRptMskFlg)
				Rpts++;
			}
		break;

	case eRPTHunmasked:		// treat any base without cRptMskFlg as a repeat masked base 
		while(Len--)
			{
			if(((Base = *pSeq++) & ~cRptMskFlg)  > eBaseT)
				return(true);
			if(!(Base & cRptMskFlg))
				Rpts++;
			}
		break;
	}

Percent = (100*Rpts) / SeqLen;
return(Percent >= MaskPercent ? true : false);
}

// LocateExacts
// Returns count of all exact matching sequences
// These matching sequences are exact (100% identity) and contain no mismatches
int	
CSeqSfx::LocateExacts(sAMProbe *pProbe,						// probe
				  unsigned int ProbeLen,					// probe length
				  etSeqBase *pProbeSeq,						// ptr to probe sequence
				  unsigned int MinHitLen,					// minimum target hit length required

  				  // following is call back on all putative hits of probe into suffix entry, allows caller to filter putative hits
				  int (* FilterCore)(sAMProbe *pProbe,				// probe
										unsigned int ProbeOfs,		// probe offset at which exact core starts 
										unsigned int ProbeCoreLen,	// length of exact core
										unsigned int TargetRef,		// user assigned target to associate with any hit to this sequence 
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]

 				  // following is call back on all accepted hits of probe into suffix entry
				  int (* ProcessCore)(sAMProbe *pProbe,			// probe
										unsigned int ProbeOfs,		// probe offset at which hit starts 
										unsigned int HitLen,		// number of bases hit
										unsigned int TargetRef,		// user assigned target to associate with any hit to this sequence
										unsigned int TargOfs,		// target loci hit	
										etSeqBase *pSeq),			// target sequence, hit starts at &pSeq[TargOfs]
				etRPTMasking RPTMasking,							// how to interpret probe bases with cRptMskFlg set
				unsigned int MaskPercent)							// if matches have more than this percent repeats then slough
{
int Rslt;

int TargPsn;	// index in pSfxArray of current suffix sequence being processed
unsigned int CurProbeOfs;
unsigned int CurMatchLen;
unsigned int CurLeftTargOfs;
unsigned int CurLeftProbeOfs;
unsigned int MinSeedLen;
unsigned int DeltaProbeOfs;

int NumProcessed;	// total number of cores after filtering
tsCPutCore *pCore;
bool bTerm;
int Idx;

tsSeqHdr *pSeqHdr;

etSeqBase *pPrb;
etSeqBase *pTrg;
etSeqBase *pPrbr;
etSeqBase *pTrgr;
etSeqBase *pMaxPrbr;
etSeqBase *pMaxTrgr;

etSeqBase *pSeq2Loc;

unsigned int CurProbeLen;

etSeqBase b1;
etSeqBase b2;

int Max2Compare;

if(MinHitLen > ProbeLen)
   return(0);

if(MinHitLen <= 12)
	{
	MinSeedLen = MinHitLen;
	DeltaProbeOfs = 1;
	}
else		// else minhitlen must be at >= 13
	{
	if(MinHitLen <= 49)
		{
		DeltaProbeOfs = MinHitLen / 4;
		MinSeedLen = MinHitLen - DeltaProbeOfs;
		if(MinSeedLen < 12)	// can only be if minhitlen < 16
			{
			MinSeedLen = 12;
			DeltaProbeOfs = MinHitLen - MinSeedLen;
			}
		}
	else	// else minhitlen must be at >= 50 so MinSeedLen will be at least 25
		{
		DeltaProbeOfs = MinHitLen / 2;
		MinSeedLen = MinHitLen - DeltaProbeOfs;
		}
	}

m_NumPutCores = 0;			// clears out any previous putative cores
NumProcessed = 0;

CurProbeOfs = 0;
bTerm = false;

do {
	pSeq2Loc = &pProbeSeq[CurProbeOfs];			// initial start
	// putative cores which contain any bases other than canonical (a,c,g or t) are simply sloughed
	if(FiltOutRepeats(eRPTHignor,0,pSeq2Loc,MinSeedLen))
		continue;

	// locate first core of at least MinSeedLen which exactly matches
	if((TargPsn = LocateFirstExact(pSeq2Loc,MinSeedLen,0,m_NumSfxEls-1)) < 0)
		continue;

	// have at least one core exactly matching for at least MinHitLen which could be extended left+right to MinHitLen
	do  {
		CurProbeLen = MinSeedLen;
		pTrg = m_pSfxEls[TargPsn].pBase;				// pts to target base
		CurLeftTargOfs = m_pSfxEls[TargPsn].SeqOfs;		// offset in target sequence

		pPrbr = &pProbeSeq[CurProbeOfs + MinSeedLen];
		pTrgr = pTrg + MinSeedLen;

		pSeqHdr = &m_pSeqHdrs[m_pSfxEls[TargPsn].SeqIdx-1];
		pMaxTrgr = pSeqHdr->pStart + pSeqHdr->SeqLen - 1; // pts to last base in probe
		pMaxPrbr = &pProbeSeq[ProbeLen-1];			// pts to last base in Seq

		// firstly try to extend left
		CurMatchLen = MinSeedLen;
		CurLeftProbeOfs = CurProbeOfs;
		if(CurLeftProbeOfs > 0 && CurLeftTargOfs > 0)
			{
			pPrb = &pProbeSeq[CurLeftProbeOfs - 1];
			pTrg -= 1;
			do	{
				b1 = (*pPrb-- & ~cRptMskFlg);
				b2 = (*pTrg-- & ~cRptMskFlg);
				if((b1 != b2) || (b1 > eBaseT || b2 > eBaseT))
					break;
				CurMatchLen++;
				CurLeftProbeOfs -= 1;
				CurLeftTargOfs -= 1;
				}
			while(CurLeftProbeOfs > 0 && CurLeftTargOfs > 0);
			}

		// after extending left then try to extend right
		while(pPrbr <= pMaxPrbr && pTrgr <= pMaxTrgr)
			{
			b1 = (*pPrbr++ & ~cRptMskFlg);
			b2 = (*pTrgr++ & ~cRptMskFlg);
			if((b1 != b2) || (b1 > eBaseT || b2 > eBaseT))
				break;
			CurMatchLen++;
			}

		if(CurMatchLen >= MinHitLen)
			{
			// apply filtering for repeats and non-canonical bases before calling AddPutCore
			if(!FiltOutRepeats(RPTMasking,MaskPercent,&pProbeSeq[CurLeftProbeOfs],CurMatchLen))
				{	
				// check if caller wants to filter this putative exact match out
				if(FilterCore != NULL)
					{
					if((Rslt = (*FilterCore)(pProbe,
								CurLeftProbeOfs,
								CurMatchLen,
								pSeqHdr->TargetRef,
								CurLeftTargOfs,
								pSeqHdr->pStart)) < 0)
						return(Rslt);
					if(!Rslt)
						{
						bTerm = true;
						continue;
						}
					}
				else
					Rslt = 1;
				if(Rslt > 0 && (Rslt =AddPutCore(CurLeftProbeOfs,	// probe offset at which hit starts
							CurMatchLen,							// number of bases hit
							pSeqHdr->SeqIdx,						// target sequence unique identifier
							CurLeftTargOfs)) < 1)					// target loci hit
					return(Rslt);
				}
			}

		// any more suffixes to check
		if(++TargPsn >= m_NumSfxEls)
			break;

		// check if probe matches next suffix
		Max2Compare = MinSeedLen;
		CurMatchLen = 0;
		pTrg = m_pSfxEls[TargPsn].pBase;
		pPrb = pSeq2Loc;
		while(Max2Compare-- && (b1 = (*pPrb++ & ~cRptMskFlg)) == (b2 = (*pTrg++ & ~cRptMskFlg)) && b1 <= eBaseT)
			CurMatchLen++;
		}
	while(!bTerm && CurMatchLen >= MinSeedLen);
	}
while(!bTerm && (CurProbeOfs += DeltaProbeOfs) <= (ProbeLen - MinSeedLen));

NumProcessed = FiltPutCores();
if(NumProcessed <= 0)
	return(NumProcessed);
pCore = m_pPutCores;
int NumAccepted = 0;
for(Idx = 0; Idx < NumProcessed; Idx++,pCore++)
	{
	if(pCore->HitLen)
		{
		if((Rslt = (*ProcessCore)(pProbe,pCore->ProbeOfs,pCore->HitLen,pSeqHdr->TargetRef,pCore->SeqOfs,pSeqHdr->pStart))<=0)
			return(Rslt);
		NumAccepted += 1;
		}
	}
return(NumAccepted);
}

int 
CSeqSfx::SortSeqSfxs( const void *arg1, const void *arg2)
{
etSeqBase *pSeq1 = ((sCSfxSeqs *)arg1)->pBase;
etSeqBase *pSeq2 = ((sCSfxSeqs *)arg2)->pBase;

etSeqBase Base1;
etSeqBase Base2;

do {
	if((Base1 = (*pSeq1++ & ~cRptMskFlg))  < (Base2 = (*pSeq2++ & ~cRptMskFlg)))
		return(-1);
	if(Base1 > Base2)
		return(1);
	}
while(Base1 == Base2 && Base1 != eBaseEOS);
return(0);
}


// SortPutCores
// Sorts putative cores by:
// ascending target header index
// ascending probe offset
// decending hit length
// ascending target seq offset
//
// BUT if hit length is 0 (indicates marked for deletion) then that core will be sorted as being last irrespective of probe and Seq offsets
int CSeqSfx::SortPutCores(const void *arg1, const void *arg2)
{
tsCPutCore *pEl1 = (tsCPutCore *)arg1;
tsCPutCore *pEl2 = (tsCPutCore *)arg2;

// hit length <= 0 are to be sorted last
if(pEl1->HitLen == 0 || pEl2->HitLen == 0)
	{
	if(pEl1->HitLen)
		return(-1);
	if(pEl2->HitLen)
		return(1);
	return(0);
	}

// sort ascending on target sequence header index
if(pEl1->SeqIdx != pEl2->SeqIdx)
	return(pEl1->SeqIdx < pEl2->SeqIdx ? -1 : 1);

// sort on ascending probe start loci
if(pEl1->ProbeOfs != pEl2->ProbeOfs)
	return(pEl1->ProbeOfs < pEl2->ProbeOfs ? -1 : 1);

// sort on decending hit length 
if(pEl1->HitLen != pEl2->HitLen)
	return(pEl1->HitLen < pEl2->HitLen ? 1 : -1);

// sort on ascending target start loci
if(pEl1->SeqOfs != pEl2->SeqOfs)
	return(pEl1->SeqOfs < pEl2->SeqOfs ? -1 : 1);

return(0);
}

