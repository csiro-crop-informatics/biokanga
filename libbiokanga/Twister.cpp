/*
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2017, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * See LICENCE for the complete licence information (https://github.com/csiro-crop-informatics/biokanga/LICENSE)
 * Contact: Alex Whan <alex.whan@csiro.au>
 */

#include "stdafx.h"
#ifdef _WIN32
#include "./commhdrs.h"
#else
#include "./commhdrs.h"
#endif



CTwister::CTwister(void)
{
m_pOctStructParams = NULL;
}

CTwister::~CTwister(void)
{
if(m_pOctStructParams!=NULL)
	delete m_pOctStructParams;
m_pOctStructParams = NULL;
}


// LoadStructParams
// Load structural parameters from file
// Basic stats are generated during the load process which covers mean,min,max and stddev for each parameter
//
int 
CTwister::LoadStructParams(char *pszStructParamsFile) // load structural parameters from file
{
int Rslt;
if((Rslt = CConformation::LoadStructOctamersParams(pszStructParamsFile)) < eBSFSuccess)
	return(Rslt);
GenStructParamStats();				// generate basic structure parameter stats
return(eBSFSuccess);
}

int 
CTwister::GenStructParamStats(void)				// generate basic structure parameter stats
{
tsOctStructParam *pStructParams;
tsStructStats *pStructStats;
int *pStructValue;
int Idx;
int ParamIdx;
INT64 Means[eSSNumStatParams];
double Variances[eSSNumStatParams];

pStructStats = m_StructParamStats;
for(ParamIdx = 0; ParamIdx < eSSNumStatParams; ParamIdx++,pStructStats++)
	{
	pStructStats->Max = INT_MIN;
	pStructStats->Min = INT_MAX;
	pStructStats->Mean= 0;
	pStructStats->StdDev=0;
	Means[ParamIdx] = 0;
	Variances[ParamIdx] = 0.0;
	}

pStructParams = m_pOctStructParams;
for(Idx = 0; Idx < cNumParamOctamers; Idx++, pStructParams++)
	{
	pStructStats = m_StructParamStats;
	for(ParamIdx = 0; ParamIdx < eSSNumStatParams; ParamIdx++,pStructStats++)
		{
		pStructValue = MapStructParam2Ptr((teOctStructStats)ParamIdx,pStructParams);
		if(*pStructValue < pStructStats->Min)
			pStructStats->Min = *pStructValue;
		if(*pStructValue > pStructStats->Max)
			pStructStats->Max = *pStructValue;
		Means[ParamIdx] += *pStructValue;
		}
	}
pStructStats = m_StructParamStats;
for(ParamIdx = 0; ParamIdx < eSSNumStatParams; ParamIdx++,pStructStats++)
	pStructStats->Mean = (int)(Means[ParamIdx]/cNumParamOctamers);

// now that means are known then stddev can be calculated
pStructParams = m_pOctStructParams;
for(Idx = 0; Idx < cNumParamOctamers; Idx++, pStructParams++)
	{
	pStructStats = m_StructParamStats;
	for(ParamIdx = 0; ParamIdx < eSSNumStatParams; ParamIdx++,pStructStats++)
		{
		pStructValue = MapStructParam2Ptr((teOctStructStats)ParamIdx,pStructParams);
		Variances[ParamIdx] += pow(((*pStructValue - pStructStats->Mean) / 100.0),2);
		}
	}

pStructStats = m_StructParamStats;
for(ParamIdx = 0; ParamIdx < eSSNumStatParams; ParamIdx++,pStructStats++)
	pStructStats->StdDev = (int)(100.0 * sqrt(Variances[ParamIdx] / cNumParamOctamers));
return(eBSFSuccess);
}

void
CTwister::SetMissing(tsOctStructParam *pRetStructParam)
{
pRetStructParam->Param.energy = INT_MIN;
pRetStructParam->Param.minorgroove = INT_MIN;
pRetStructParam->Param.majorgroove = INT_MIN;
pRetStructParam->Param.rise = INT_MIN;
pRetStructParam->Param.rmsd = INT_MIN;
pRetStructParam->Param.roll = INT_MIN;
pRetStructParam->Param.shift = INT_MIN;
pRetStructParam->Param.slide = INT_MIN;
pRetStructParam->Param.tilt = INT_MIN;
pRetStructParam->Param.twist = INT_MIN;
pRetStructParam->Param.orchid = INT_MIN;
}


// GetStructParams
// Get structural parameters for specified step in sequence
// Steps are numbered 1..SeqLen-1, with step 1 between 1st and 2nd base of sequence
// SeqLen must be a minimum of 8
// Structural parameters will be set to INT_MIN if any nucleotide in the octamer is not one of a,c,g,t
// 
int
CTwister::GetStructParams(unsigned int Step,	// which step (1..SeqLen-1)
				unsigned int SeqLen,			// sequence length (8..n)
				etSeqBase *pSeq,				// sequence
				tsOctStructParam *pRetStructParam) // where to return structural parameters
{
int OctIdx;
etSeqBase *pOctamer;
etSeqBase octamer[8];
int	Cnt;
int	OctOfs;
int	SeqOfs;

if(SeqLen < 8 || Step < 1 || Step >= SeqLen || pSeq == NULL || pRetStructParam == NULL || m_pOctStructParams == NULL)
	return(eBSFerrParams);

// if can use the midstep of an octamer then do so...
if(Step > 3 && Step < SeqLen - 3)
	{
	pOctamer = pSeq + Step - 4;
	if((OctIdx = StructParamIdx(pOctamer)) < 0)
		SetMissing(pRetStructParam);
	else
		*pRetStructParam = m_pOctStructParams[OctIdx];
	return(eBSFSuccess);
	}

// not a octamer midstep, need to interpolate
// step    cnt   toofs fromofs
// 1	   5     3		3
// 2       6     2		2
// 3       7     1      1
//
// 5       7     0		SeqLen - 8
// 6       6	 0      SeqLen - 7
// 7       5	 0		SeqLen - 6
memset(octamer,eBaseA,8);
if(Step < 4)
	{
	Cnt = 4 + Step;
	OctOfs = 4 - Step;
	SeqOfs = 0;
	}
else
	{
	Step = 8 - (SeqLen - Step);						// adjust Step to 5..7
	Cnt = 12 - Step;
	OctOfs = 0;
	SeqOfs = SeqLen - (Cnt + 1);
	}
memmove(&octamer[OctOfs],&pSeq[SeqOfs],Cnt);

Interpolate(Step,octamer,pRetStructParam);
return(eBSFSuccess);
}

int
CTwister::GetStructParam(unsigned int Step,	// which step (1..SeqLen-1)
				unsigned int SeqLen,		// sequence length (8..n)
				etSeqBase *pSeq,			// sequence
				teOctStructStats Param)		// which structural parameter value to return
{
tsOctStructParam StructParams;
int *pStructValue;
if(GetStructParams(Step,SeqLen,pSeq,&StructParams)!=eBSFSuccess)
	{
	AddErrMsg("CTwister::GetStructParam","GetStructParams failed for Step: %d SeqLen: %d",Step,SeqLen);
	return(INT_MIN);
	}
pStructValue = MapStructParam2Ptr(Param,&StructParams);
return(*pStructValue);
}

int
CTwister::Interpolate(unsigned int Step,	// which step to interpolate (1..7)
			etSeqBase *pSeq,				// known octamer sequence left (step 5..7) or right (step 1..3) filled with eBaseA
			tsOctStructParam *pRetStructParam) // where to returned interpolated structural parameters
{
int energy,minorgroove,majorgroove,twist,roll,tilt,rise,slide,shift,rmsd;
tsOctStructParam *pStruct;
int Cnt;
int Iters;

etSeqBase octamer[8];						
int OctIdx;
bool bLeft;
if(Step < 1 || Step > 7 || pSeq == NULL || pRetStructParam == NULL || m_pOctStructParams == NULL)
	return(eBSFerrParams);

if(Step == 4)								// if octamer midstep then no need to interpolate!
	return(GetStructParams(4,8,pSeq,pRetStructParam));

memmove(octamer,pSeq,8);						// local copy because we will be modifying it..

if(Step > 4)								// sequence to left is known
	{
	bLeft = true;
	Step -= 4;						
	}
else										// else sequence to right is known
	{
	bLeft = false;
	Step = 4 - Step;	
	}

switch(Step) {
	case 1: 
		Iters = 4;
		break;
	case 2:
		Iters = 16;
		break;
	default:
		Iters = 64;
		break;
	}

energy=minorgroove=majorgroove=twist=roll=tilt=rise=slide=shift=rmsd=0;

for(Cnt = 0; Cnt < Iters; Cnt++)
	{
	OctIdx = StructParamIdx(octamer);
	if(OctIdx < 0)
		{
		SetMissing(pRetStructParam);
		return(eBSFSuccess);
		}

	pStruct = &m_pOctStructParams[OctIdx];
	if(pStruct->Param.energy == 0)
		return(eBSFerrStructParm);

	energy += pStruct->Param.energy;
	minorgroove += pStruct->Param.minorgroove;
	majorgroove += pStruct->Param.majorgroove;
	twist += pStruct->Param.twist;
	roll += pStruct->Param.roll;
	tilt += pStruct->Param.tilt;
	rise += pStruct->Param.rise;
	slide += pStruct->Param.slide;
	shift += pStruct->Param.shift;
	rmsd += pStruct->Param.rmsd;
	
	if(!bLeft)
		{
		octamer[0]++;
		if(octamer[0] > eBaseT)
			{
			octamer[0] = eBaseA;
			octamer[1]++;
			if(octamer[1] > eBaseT)
				{
				octamer[1] = eBaseA;
				octamer[2]++;
				}
			}
		}
	else
		{
		octamer[7]++;
		if(octamer[7] > eBaseT)
			{
			octamer[7] = eBaseA;
			octamer[6]++;
			if(octamer[6] > eBaseT)
				{
				octamer[6] = eBaseA;
				octamer[5]++;
				}
			}
		}
	}
pRetStructParam->Param.energy = energy / Iters;
pRetStructParam->Param.minorgroove = minorgroove / Iters;
pRetStructParam->Param.majorgroove = majorgroove / Iters;
pRetStructParam->Param.twist = twist / Iters;
pRetStructParam->Param.roll = roll / Iters;
pRetStructParam->Param.tilt = tilt / Iters;
pRetStructParam->Param.rise = rise / Iters;
pRetStructParam->Param.slide = slide / Iters;
pRetStructParam->Param.shift = shift / Iters;
pRetStructParam->Param.rmsd = rmsd / Iters;
return(eBSFSuccess);
}

// determine idex to use for assumed octamer sequence
// eBSFerrStructStep returned if any base in octamer is indeterminate - 'N'
int
CTwister::StructParamIdx(etSeqBase *pOctamer)		// sequence
{
etSeqBase Base;
int OctIdx = 0;
int Len = 8;
if(pOctamer == NULL || m_pOctStructParams == NULL)
	return(eBSFerrParams);

while(Len--)
	{
	OctIdx <<= 2;
	Base = *pOctamer++  & ~cRptMskFlg;
	if(Base > eBaseT)
		return(eBSFerrStructStep);								// unrecognised base
	OctIdx |= Base;
	}

return(OctIdx);
}

int
CTwister::ProcessSequence(int hRslts,			// file to write sequence structure into
				  unsigned int RefID,			// non-unique identifier by which related sequences can be later associated 
				  char *pszSpecies,				// species description
				  unsigned int EntryID,			// sequence entry identifier (could be chromosome) for this sequence
				  char *pszEntry,				// entry or chromosome description
				  unsigned int EntryOfs,		// entry or chromosome offset (0..n) corresponding to left base of 1st step
				  unsigned int iStartOfs,		// initial starting offset (0..n) in pSeq
				  unsigned int iNumSteps,		// number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  unsigned int SeqLen,			// total length of sequence
				  etSeqBase *pSeq,				// sequence to be processed
				  bool bXML)					// results file type: false == CSV, true == XML
{
char szVal[80];
char szLineBuff[2000];
int BuffLen;
tsOctStructParam StructParam;
unsigned int Step;
unsigned int LastStep;
int iNumDP;
int Rslt;

if(!EntryID ||  SeqLen < 8 || iStartOfs >= SeqLen - 1 || 
   iStartOfs + iNumSteps >= SeqLen ||
   pSeq == NULL || m_pOctStructParams == NULL)
	return(eBSFerrParams);

if(hRslts < 0)
	{
	AddErrMsg("CTwister::ProcessSequence","No file opened to write results into");
	return(eBSFerrClosed);
	}

if(iNumSteps == 0)
	iNumSteps = SeqLen - iStartOfs - 1;

LastStep = iStartOfs + iNumSteps;
iNumDP = 0;
Rslt = eBSFSuccess;
for(Step = iStartOfs+1; Step <= LastStep && Rslt == eBSFSuccess; Step++)
	{
	if((Rslt=GetStructParams(Step,SeqLen,pSeq,&StructParam))==eBSFSuccess)
		{
		if(!bXML)
			{
			BuffLen = sprintf(szLineBuff,"\"%d\",",RefID);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\"%s\",",pszSpecies);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\"%s\",",pszEntry);
			BuffLen += sprintf(&szLineBuff[BuffLen],"%d,",EntryID);
			BuffLen += sprintf(&szLineBuff[BuffLen],"%d,",EntryOfs++);
			BuffLen += sprintf(&szLineBuff[BuffLen],"%d,",Step);
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.minorgroove,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.majorgroove,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.energy,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.twist,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.roll,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.tilt,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.rise,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.slide,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.shift,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s,",Fmt2FixedDec(StructParam.Param.rmsd,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"%s\n",Fmt2FixedDec(StructParam.Param.orchid,szVal));
			}
		else
			{
					// output results here as XML
			BuffLen = sprintf(szLineBuff,"<Step>\n");
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<refid>%d</refid>\n",RefID);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<species>%s</species>\n",pszSpecies);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<entry>%s</entry>\n",pszEntry);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<entryid>%d</entryid>\n",EntryID);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<entryofs>%d</entryofs>\n",EntryOfs++);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<step>%d</step>\n",Step);
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<minor>%s</minor>\n",Fmt2FixedDec(StructParam.Param.minorgroove,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<major>%s</minor>\n",Fmt2FixedDec(StructParam.Param.majorgroove,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<energy>%s</energy>\n",Fmt2FixedDec(StructParam.Param.energy,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<twist>%s</twist>\n",Fmt2FixedDec(StructParam.Param.twist,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<roll>%s</roll>\n",Fmt2FixedDec(StructParam.Param.roll,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<tilt>%s</tilt>\n",Fmt2FixedDec(StructParam.Param.tilt,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<rise>%s</rise>\n",Fmt2FixedDec(StructParam.Param.rise,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<slide>%s</slide>\n",Fmt2FixedDec(StructParam.Param.slide,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<shift>%s</shift>\n",Fmt2FixedDec(StructParam.Param.shift,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<rmsd>%s</rmsd>\n",Fmt2FixedDec(StructParam.Param.rmsd,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"\t<orchid>%s</orchid>\n",Fmt2FixedDec(StructParam.Param.orchid,szVal));
			BuffLen += sprintf(&szLineBuff[BuffLen],"</Step>\n");
			}
		if(write(hRslts,szLineBuff,BuffLen)!=BuffLen)
			{
			AddErrMsg("CTwister::ProcessSequence","Unable to write to results file: %s\n",strerror(errno));
			return(eBSFerrFileAccess);
			}
		}
	}
return(Rslt);
}

int
CTwister::GetSequenceConformation(teOctStructStats confparam,	// process for this conformational parameter
				  int iStartOfs,			// initial starting offset (0..n) in pSeq
				  int iNumSteps,			// number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  int SeqLen,				// total length of sequence
				  etSeqBase *pSeq,			// sequence to be processed
				  int *pRetValues)			// where to return step conformational values
{
tsOctStructParam StructParam;
unsigned int Step;
unsigned int LastStep;
int iNumDP;
int Rslt;

if(SeqLen < 8 || iStartOfs >= SeqLen - 1 || 
   iStartOfs + iNumSteps >= SeqLen || pRetValues == NULL ||
   pSeq == NULL || m_pOctStructParams == NULL || confparam >= eSSNumStatParams)
	return(eBSFerrParams);

if(iNumSteps == 0)
	iNumSteps = SeqLen - iStartOfs - 1;

LastStep = iStartOfs + iNumSteps;
iNumDP = 0;
Rslt = eBSFSuccess;
for(Step = iStartOfs+1; Step <= LastStep && Rslt == eBSFSuccess; Step++)
	{
	if((Rslt=GetStructParams(Step,SeqLen,pSeq,&StructParam))==eBSFSuccess)
		{
		switch(confparam) {
			 case eSSenergy:		// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
				 *pRetValues++ = StructParam.Param.energy;
				 break;

			 case eSSminorgroove:				// minor groove int (dimensions * 10000) e.g 10.784 ==> 107840
				 *pRetValues++ = StructParam.Param.minorgroove;
				 break;

			 case eSSmajorgroove:				// major groove int (dimensions * 10000) e.g 10.784 ==> 107840
				 *pRetValues++ = StructParam.Param.majorgroove;
				 break;

			 case eSStwist:					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
				 *pRetValues++ = StructParam.Param.twist;
				 break;

			 case eSSroll:					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
				 *pRetValues++ = StructParam.Param.roll;
				 break;

			 case eSStilt:					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
				 *pRetValues++ = StructParam.Param.tilt;
				 break;

			 case eSSrise:					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
				 *pRetValues++ = StructParam.Param.rise;
				 break;

			 case eSSslide:					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
				 *pRetValues++ = StructParam.Param.slide;
				 break;

			 case eSSshift:					// shift int(shift * 10000) e.g 	0.0645 ==> 645
				 *pRetValues++ = StructParam.Param.shift;
				 break;

			 case eSSrmsd:					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
				 *pRetValues++ = StructParam.Param.rmsd;
				 break;

			 case eSSORChidVal:				// hydroxyl radical cleavage value from ORChid dataset
				 *pRetValues++ = StructParam.Param.orchid;
				 break;
			}
		}
	}
return(Rslt);
}


char *
CTwister::Fmt2FixedDec(int Value,char *pszRet)
{
static char szFmtBuff[8];
int Intg = Value/10000;
int Remainder = abs(Value)%10000;
if(pszRet == NULL)
	pszRet = szFmtBuff;
sprintf(pszRet,"%d.%04d",Intg,Remainder);
return(pszRet);
}



// CalcTwistStats
// Calculates the mean, min, max, and std deviation for all dinucleotide base steps in specified sequence
int	
CTwister::CalcTwistStats(unsigned int iStartOfs,		// initial starting offset (0..n) in pSeq
					unsigned int iNumSteps,		// number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				    unsigned int SeqLen,		// total length of sequence
					etSeqBase *pSeq,			// start of sequence
					teOctStructStats StructParam,	// which structural parameter to gen stats for
					tsStructStats *pStats)		// returned stats
{
tsOctStructParam StructParams;
int *pStructValue;
INT64 SumValue;
int MeanValue;
int StdDev;
unsigned int Step;
unsigned int LastStep;

if( SeqLen < 8 || iStartOfs >= SeqLen - 1 || 
   iStartOfs + iNumSteps >= SeqLen ||
   pSeq == NULL ||
   pStats == NULL)
	{
	AddErrMsg("CTwister::CalcTwistStats","One or more invalid parameters passed to function");
	return(eBSFerrParams);
	}

pStructValue = MapStructParam2Ptr(StructParam,&StructParams);
if(pStructValue == NULL)
	{
	AddErrMsg("CTwister::CalcTwistStats","Unrecognised structural parameter %d",StructParam);
	return(eBSFerrStructParm);
	};

if(iNumSteps == 0)
	iNumSteps = SeqLen - iStartOfs - 1;
LastStep = iStartOfs + iNumSteps;

SumValue = 0;

int ValidSteps = 0;
for(Step = iStartOfs+1; Step <= LastStep; Step++)
	{
	GetStructParams(Step,SeqLen,pSeq,&StructParams);
	if(StructParams.Param.energy != INT_MIN)
		{
		ValidSteps++;
		SumValue += *pStructValue;
		}
	}
if(!ValidSteps)
	{
	pStats->Mean = 0;
	pStats->Max  = 0;
	pStats->Min  = 0;
	pStats->StdDev = 0;
	return(eBSFerrStructStep);
	}

MeanValue = (int)(SumValue/(INT64)(ValidSteps));

// now that mean is known then can calculate stddev
double Variance = 0.0;
pStats->Max = INT_MIN;
pStats->Min = INT_MAX;
ValidSteps = 0;
for(Step = iStartOfs+1; Step <= LastStep; Step++)
	{
	GetStructParams(Step,SeqLen,pSeq,&StructParams);
	if(StructParams.Param.energy != INT_MIN)
		{
		ValidSteps++;
		if(*pStructValue < pStats->Min)
			pStats->Min = *pStructValue;
		if(*pStructValue > pStats->Max)
			pStats->Max = *pStructValue;
		Variance += pow(((*pStructValue - MeanValue) / 100.0),2);
		}
	}

StdDev = (int)(100.0 * sqrt(Variance / ValidSteps));
pStats->Mean = MeanValue;
pStats->StdDev= StdDev;
	
return(eBSFSuccess);
}


int
CTwister::CalcStats(char *pszSpecies,			// species - e.g 'hg17'
				  unsigned int EntryID,			// sequence entry identifier for this sequence
				  char *pszChrom,				// from this chromosome
				  unsigned int iStartOfs,		// initial starting offset (0..n) in pSeq 
				  unsigned int iNumSteps,		// number of steps to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				  unsigned int SeqLen,			// total length of sequence
				  etSeqBase *pSeq,				// sequence
				  teOctStructStats StructParam,	// which structural parameter to gen stats for
  				  unsigned int SlideWinSteps, 	// sliding window size in dinucleotide steps (5..100000)
				  int hRslts,					// results file to write into
					bool bXML)					// results file type: false == CSV, true == XML
{
tsStructStats HistBins[cNumHistBins];			// to hold binned step structural value histogram
int BinIdx;
unsigned int LastStep, Step, WinSteps;
tsStructStats CurStructStats;						// stats for current window

if(pszSpecies == NULL || pszSpecies[0] == '\0' ||
   pszChrom == NULL || pszChrom[0] == '\0' ||
   !EntryID || SeqLen < 8 || iStartOfs >= SeqLen - 1 || iStartOfs + iNumSteps >= SeqLen ||
   pSeq == NULL || hRslts < 0)
	{
	AddErrMsg("CTwister::CalcStats","One or more invalid parameters passed to function");
	return(eBSFerrParams);
	}



// limit WinNumSteps to be between 1 and 100000 inclusive
if(iNumSteps < SlideWinSteps)
	SlideWinSteps = iNumSteps;

if(SlideWinSteps > 100000)
	SlideWinSteps = 100000;
else
	if(SlideWinSteps < 1)
		SlideWinSteps = 1;

memset(HistBins,0,sizeof(HistBins));

if(iNumSteps == 0)
	iNumSteps = SeqLen - iStartOfs - 1;
LastStep = iStartOfs + iNumSteps;
WinSteps = SlideWinSteps;
for(Step = iStartOfs; Step < LastStep; Step++)
	{
	if(SeqLen - Step <= WinSteps)
		WinSteps = SeqLen - Step - 1;
	CalcTwistStats(Step,WinSteps,SeqLen,pSeq,StructParam,&CurStructStats);
	HistBins[CurStructStats.Mean].Mean++;
	HistBins[CurStructStats.StdDev].StdDev++;
	HistBins[CurStructStats.Max].Max++;
	HistBins[CurStructStats.Min].Min++;
	}

	// histogram has been calculated...
	// output bins from BinMinIdx to BinMaxIdx
for(BinIdx = 0; BinIdx < cNumHistBins; BinIdx++)
	{
	if(HistBins[BinIdx].Min > 0)
		OutResult(bXML,hRslts,pszSpecies,pszChrom,(char *)"Min",EntryID,0,SlideWinSteps,BinIdx,HistBins[BinIdx].Min);
	if(HistBins[BinIdx].Max > 0)
		OutResult(bXML,hRslts,pszSpecies,pszChrom,(char *)"Max",EntryID,1,SlideWinSteps,BinIdx,HistBins[BinIdx].Max);
	if(HistBins[BinIdx].Mean > 0)
		OutResult(bXML,hRslts,pszSpecies,pszChrom,(char *)"Mean",EntryID,2,SlideWinSteps,BinIdx,HistBins[BinIdx].Mean);
	if(HistBins[BinIdx].StdDev > 0)
		OutResult(bXML,hRslts,pszSpecies,pszChrom,(char *)"StdDev",EntryID,3,SlideWinSteps,BinIdx,HistBins[BinIdx].StdDev);
	}


return(eBSFSuccess);
}


void
CTwister::OutResult(bool bXML,				// true==XML, false==CSV
					int hRslts,				// file to output into
					char *pszSpecies,		// species name e.g hs17
					char *pszChrom,			// chromosome e.g chr1
					char *pszStat,			// statistics type e.g mean
					int EntryID,			// entry ident corresponding to pszChrom
					int Stat,				// statistics type: 0 == Min, 1== Max, 2==Mean, 4==StdDev
					int Steps,				// number of dinucleotide steps stats generated over (1..n)
					int Value,				// stats value
					int Freq)				// frequency or number of instances of Value
{
char szBuff[2000];
int BuffL;
if(!bXML)
	BuffL = sprintf(szBuff,"'%s','%s','%s',%d,%d,%d,%d,%d\n",pszSpecies,pszChrom,pszStat,EntryID,Stat,Steps,Value,Freq);
else
	{
	BuffL = sprintf(szBuff,"<bin>\n\t<species>%s</species><chrom>%s</chrom><stat>%s</stat>",pszSpecies,pszChrom,pszStat);
	BuffL += sprintf(&szBuff[BuffL],"\n\t<chromid>%d</chromid><statid>%d</statid><steps>%d</steps><value>%d</value><freq>%d</freq>\n</bin>\n",
				EntryID,Stat,Steps,Value,Freq);
	}
CUtility::SafeWrite(hRslts,szBuff,BuffL);
}

int
CTwister::CalcStruct(char *pszBioSeqFile,	// sequence file
				char *pszSpecies,			// reference species name
				int iSeqEntry,				// which entry (1..n) or all entries (0) in sequence file
				int iStartOfs,				// process from start position (0..n)
				int  iNumSteps,				// process this many steps starting at Seq[StartOfs]|Seq[StartOfs+1]
				char *pszResultsFile,		// where to write results
				bool bXML,					// results file type: false == CSV, true == XML
				char *pszStructParamsFile)	// structural parameters file, NULL if to use existing loaded
{
unsigned int ProbeEntryID;
CBioSeqFile *pBioSeqFile;
int hRslts = -1;
unsigned int AllocLen = 0;
unsigned int SeqLen;
unsigned int iRemaining;
unsigned char *pSeqBuff = NULL;
char szRlstsFile[_MAX_PATH];
unsigned char szSource[cBSFSourceSize];
unsigned char szDescription[cBSFDescriptionSize];
char szProbeDescr[cBSFSourceSize+cBSFDescriptionSize+1];
int LBuff;
char Buff[100];
int Rslt;


strcpy(szRlstsFile,pszResultsFile);
if(pszStructParamsFile != NULL && pszStructParamsFile[0] != '\0')
	{
	if((Rslt=LoadStructParams(pszStructParamsFile))!=eBSFSuccess)
		{
		AddErrMsg("CTwister::CalcStruct","Unable to load structure parameters from %s",pszStructParamsFile);
		return(Rslt);
		}
	}
else
	{
	if(m_pOctStructParams == NULL)
		{
		AddErrMsg("CTwister::CalcStruct","Structural parameter file not specified but no structural parameters previously loaded");
		return(eBSFerrStructParm);
		}
	}

if((pBioSeqFile = new CBioSeqFile) == NULL)
	{
	AddErrMsg("CTwister::CalcStruct","Unable to instantiate CBioSeqFile");
	return(eBSFerrMem);
	}

if((Rslt=pBioSeqFile->Open(pszBioSeqFile,cBSFTypeSeq,false))!=eBSFSuccess)
	{
	delete pBioSeqFile;
	return(Rslt);
	}

#ifdef _WIN32
hRslts = open(pszResultsFile, O_CREATETRUNC);
#else
if((hRslts = open(pszResultsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!= -1)
  if(ftruncate(hRslts,0)!=0)
		{
		AddErrMsg("CTwister::CalcStruct","Unable to truncate %s - %s",pszResultsFile,strerror(errno));
		Rslt =eBSFerrCreateFile;
		return(Rslt);
		}
#endif
if(hRslts < 0)
	{
	AddErrMsg("CTwister::CalcStruct","Unable to create or truncate results file %s error: %s",pszResultsFile,strerror(errno));
	delete pBioSeqFile;
	return(eBSFerrCreateFile);
	}
if(bXML) 
	{
	LBuff = sprintf(Buff,"<dataset>\n");
	if(write(hRslts,Buff,LBuff) != LBuff)
		{
		AddErrMsg("CTwister::CalcStruct","Unable to write to results file %s error: %s",pszResultsFile,strerror(errno));
		close(hRslts);
		delete pBioSeqFile;
		return(eBSFerrFileAccess);
		}
	}

if(iSeqEntry)
	{
	if((Rslt=pBioSeqFile->Exists(iSeqEntry))!=eBSFSuccess)
		{
		AddErrMsg("CTwister::CalcStruct","Unable to locate specified entry %d in probe %s",iSeqEntry,pszBioSeqFile);
		close(hRslts);
		delete pBioSeqFile;
		return(Rslt);
		}
	}



if(iNumSteps == 0)
	iNumSteps = (unsigned int )-1;

ProbeEntryID = iSeqEntry;
while(iSeqEntry || (ProbeEntryID = pBioSeqFile->Next(ProbeEntryID)) > eBSFSuccess)
		{
		pBioSeqFile->GetNameDescription(ProbeEntryID,cBSFSourceSize-1,(char *)&szSource,
								cBSFDescriptionSize-1,(char *)&szDescription);
		sprintf(szProbeDescr,"%s|%s",szSource,szDescription);

		SeqLen = pBioSeqFile->GetDataLen(ProbeEntryID);
		if(!SeqLen)
			continue;

		if((unsigned int)iStartOfs >= SeqLen)
			continue;

		if(AllocLen < (SeqLen + 10))
			{
			if(pSeqBuff != NULL)
				delete pSeqBuff;
			AllocLen = 0;
			if((pSeqBuff = new unsigned char [SeqLen + 10]) == NULL)
				{
				AddErrMsg("CTwister::CalcStruct","Unable to locate alloc memory for pSeqBuff %d bytes",SeqLen + 10);
				close(hRslts);
				return(eBSFerrMem);
				}
			AllocLen = SeqLen + 10;
			}
		if(pBioSeqFile->GetData(ProbeEntryID,eSeqBaseType,0,pSeqBuff,SeqLen) != SeqLen)
			break;
		pSeqBuff[SeqLen] = eBaseEOS;	// right extensions depend on sequence being terminated!

		iRemaining = SeqLen - iStartOfs;
		if(iRemaining > (unsigned)iNumSteps)
			iRemaining = iNumSteps;

		ProcessSequence(hRslts,								// file to write sequence structure into
						ProbeEntryID,
						pszSpecies,						// species
						ProbeEntryID,					// sequence entry identifier for this sequence
						szProbeDescr,					// from this chromosome
						 iStartOfs,						// chomosome psn
						  iStartOfs,						// where in sequence to start (0..SeqLen - 8)
						  iRemaining,						// number of mid octamer steps to process
						  SeqLen,							// total length of sequence
						  pSeqBuff,
			  			  bXML);
#ifdef _WIN32
		_commit(hRslts);
#else
		fsync(hRslts);
#endif
		if(iSeqEntry)
			break;
		}

if(pSeqBuff != NULL)
	delete pSeqBuff;

if(hRslts > 0)
	{
	if(bXML)
		{
		LBuff = sprintf(Buff,"</dataset>\n");
		if(write(hRslts,Buff,LBuff) != LBuff)
			{
			AddErrMsg("CTwister::CalcStruct","Unable to write to results file %s error: %s",pszResultsFile,strerror(errno));
			delete pBioSeqFile;
			close(hRslts);
			return(eBSFerrFileAccess);
			}
		}
	close(hRslts);
	}

delete pBioSeqFile;
return(eBSFSuccess);
}


// GenRefTwistDist
int
CTwister::GenRefStructDist(char *pszResultsFile,				  // file to contain results
				bool bXML,					// results file type: false == CSV, true == XML
				char *pszStructParamsFile) // structural parameters file, NULL to use existing loaded file
{
char szBuff[2000];
int BuffL;
char szVal[80];
tsOctStructParam *pStructParams;
tsStructStats *pStructStats;
int *pStructValue;
int hRslts;
int LBuff;
int BinRange[eSSNumStatParams];
int HistBinIdx;

int HistBins[eSSNumStatParams][cNumHistBins];
CBioSeqFile BioSeqFile;
unsigned int Idx,ParamIdx;
char *pszParam;
int Rslt;

if(pszStructParamsFile != NULL && pszStructParamsFile[0] != '\0')
	{
	if((Rslt=LoadStructParams(pszStructParamsFile))!=eBSFSuccess)
		{
		AddErrMsg("CTwister::GenRefStructDist","Unable to load structure parameters from %s",pszStructParamsFile);
		return(Rslt);
		}
	}
else
	{
	if(m_pOctStructParams == NULL)
		{
		AddErrMsg("CTwister::GenRefStructDist","Structural parameter file not specified but no structural parameters previously loaded");
		return(eBSFerrStructParm);
		}
	}

#ifdef _WIN32
hRslts = open(pszResultsFile,O_CREATETRUNC);
#else
if((hRslts = open(pszResultsFile,O_RDWR | O_CREAT,S_IREAD | S_IWRITE))!=-1)
	  if(ftruncate(hRslts,0)!=0)
			{
			AddErrMsg("CTwister::GenRefStructDist","Unable to truncate %s - %s",pszResultsFile,strerror(errno));
			Rslt =eBSFerrCreateFile;
			return(Rslt);
			}
#endif
if(hRslts < 0)
		{
		AddErrMsg("CTwister::GenRefStructDist","Unable to create or truncate results file %s error: %s",pszResultsFile,strerror(errno));
		return(eBSFerrCreateFile);
		}

if(bXML)
	{
	LBuff = sprintf(szBuff,"<dataset>\n");
	if(write(hRslts,szBuff,LBuff) != LBuff)
		{
			AddErrMsg("CTwister::GenRefStructDist","Unable to write to results file %s error: %s",pszResultsFile,strerror(errno));
			close(hRslts);
			return(eBSFerrFileAccess);
		}
	}

pStructStats = m_StructParamStats;
for(ParamIdx = 0; ParamIdx < eSSNumStatParams; ParamIdx++,pStructStats++)
	if(pStructStats->Max == pStructStats->Min)
		BinRange[ParamIdx] = 0;
	else
		BinRange[ParamIdx] = (pStructStats->Max - pStructStats->Min - 1)/cNumHistBins;

memset(HistBins,0,sizeof(HistBins));
pStructParams = m_pOctStructParams;
for(Idx = 0; Idx < cNumParamOctamers; Idx++, pStructParams++)
	{
	pStructStats = m_StructParamStats;
	for(ParamIdx = 0; ParamIdx < eSSNumStatParams - 1; ParamIdx++,pStructStats++)
		{
		pStructValue = MapStructParam2Ptr((teOctStructStats)ParamIdx,pStructParams);
		if(BinRange[ParamIdx] != 0)
			HistBinIdx = (*pStructValue - pStructStats->Min)/BinRange[ParamIdx];
		else
			HistBinIdx = 0;
		HistBins[ParamIdx][HistBinIdx] +=1;
		}
	}
		// histogram has been calculated...
		// output bins from BinMinIdx to BinMaxIdx
pStructStats = m_StructParamStats;
for(ParamIdx = 0; ParamIdx < eSSNumStatParams-1; ParamIdx++,pStructStats++)
	{
	pszParam = (char *)MapStructParam2Txt((teOctStructStats)ParamIdx);
	for(HistBinIdx = 0; HistBinIdx < cNumHistBins; HistBinIdx++)
		{
		if(HistBins[ParamIdx][HistBinIdx] > 0)
			{
			if(!bXML)
				BuffL = sprintf(szBuff,"%s,%d,%d,%s,%d\n",pszParam,ParamIdx,HistBinIdx,
							Fmt2FixedDec((HistBinIdx * BinRange[ParamIdx]) + pStructStats->Min,szVal),
							HistBins[ParamIdx][HistBinIdx]);
			else
				{
				BuffL = sprintf(szBuff,"<histbin>\n\t<param>%s</param><paramid>%d</paramid><bin>%d</bin>",
													pszParam,ParamIdx,HistBinIdx);
				BuffL += sprintf(&szBuff[BuffL],"<value>%s</value><freq>%d</freq>\n</histbin>\n",
					Fmt2FixedDec((HistBinIdx * BinRange[ParamIdx]) + pStructStats->Min,szVal),
							HistBins[ParamIdx][HistBinIdx]);

				}
			if(write(hRslts,szBuff,BuffL)!= BuffL)
				{
				AddErrMsg("CTwister::GenRefStructDist","Unable to write to results file %s error: %s",pszResultsFile,strerror(errno));
				close(hRslts);
				return(eBSFerrFileAccess);
				}
			}
		}
	}


if(bXML && hRslts != -1)
	{
	LBuff = sprintf(szBuff,"</dataset>\n");
	if(write(hRslts,szBuff,LBuff) != LBuff)
		{
		AddErrMsg("CTwister::GenRefStructDist","Unable to write to results file %s error: %s",pszResultsFile,strerror(errno));
		close(hRslts);
		return(eBSFerrFileAccess);
		}
	}
if(hRslts != -1)
	close(hRslts);
return(eBSFSuccess);
}


int	
CTwister::CalcDiffStats(unsigned int iRefStartOfs,	// initial starting offset (0..n) in pRefSeq
					unsigned int iNumSteps,		// number of steps to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
				    unsigned int RefSeqLen,		// total length of reference sequence
					etSeqBase *pRefSeq,			// start of reference sequence
					unsigned int iRelStartOfs,	// initial starting offset (0..n) in pRelSeq
					unsigned int RelSeqLen,		// total length of reference sequence
					etSeqBase *pRelSeq,			// start of relative sequence
					teOctStructStats StructParam,	// which structural parameter to gen stats for
					tsStructStats *pStats)		// returned stats
{
unsigned int Step;
unsigned int LastStep;
int RefValue, RelValue;
int StdDev, Min, Max;
INT64 Mean;
int Diff;
double Variance;

if(iNumSteps < 7)
	{
	AddErrMsg("CTwister::CalcDiffStats","Number of steps must be 7 or more");
	return(eBSFerrParams);
	}

if(RefSeqLen < 8 || iRefStartOfs >= RefSeqLen - 1 || 
   iRefStartOfs + iNumSteps >= RefSeqLen ||
   pRefSeq == NULL || pStats == NULL)
	{
	AddErrMsg("CTwister::CalcDiffStats","One or more invalid ref parameters passed to function");
	return(eBSFerrParams);
	}

if(RelSeqLen < 8 || iRelStartOfs >= RelSeqLen - 1 || 
   iRelStartOfs + iNumSteps >= RelSeqLen ||
   pRelSeq == NULL)
	{
	AddErrMsg("CTwister::CalcDiffStats","One or more invalid rel parameters passed to function");
	return(eBSFerrParams);
	}
Variance = 0.0;
LastStep = iNumSteps;
Max = INT_MIN;
Min = INT_MAX;
Mean= 0;

for(Step = 1; Step <= LastStep; Step++)
	{
	RefValue = GetStructParam(Step + iRefStartOfs,RefSeqLen,pRefSeq,StructParam);
	RelValue = GetStructParam(Step + iRelStartOfs,RefSeqLen,pRelSeq,StructParam);
	Diff = RefValue - RelValue;
	Variance += pow(((Diff) / 100.0),2);
	if(pStats != NULL)
		{
		Diff = abs(Diff);
		if(Diff > Max)
			Max = Diff;
		if(Diff < Min)
			Min = Diff;
		Mean += Diff;
		}
	}

StdDev = (int)(100.0 * sqrt(Variance / iNumSteps));

pStats->Max = Max;
pStats->Min = Min;
pStats->Mean= (int)(Mean / iNumSteps);
pStats->StdDev = StdDev;

return(eBSFSuccess);
}


int *
CTwister::MapStructParam2Ptr(teOctStructStats Param,	// which parameter
							 tsOctStructParam *pStruct) // which parameters instance
{
switch(Param) {
	case eSSenergy:				// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
		return(&pStruct->Param.energy);

	case eSSminorgroove:				// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		return(&pStruct->Param.minorgroove);

	case eSSmajorgroove:				// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		return(&pStruct->Param.majorgroove);

	case eSStwist:					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
		return(&pStruct->Param.twist);

	case eSSroll:					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
		return(&pStruct->Param.roll);

	case eSStilt:					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
		return(&pStruct->Param.tilt);

	case eSSrise:					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
		return(&pStruct->Param.rise);

	case eSSslide:					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
		return(&pStruct->Param.slide);

	case eSSshift:					// shift int(shift * 10000) e.g 	0.0645 ==> 645
		return(&pStruct->Param.shift);

	case eSSrmsd:					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
		return(&pStruct->Param.rmsd);

	case eSSORChidVal:
		return(&pStruct->Param.orchid);

	default:
		return(NULL);
	};
}

const char *
CTwister::MapStructParam2Txt(teOctStructStats Param)
{
const char *pszParam;
switch(Param) {
	case eSSenergy:				// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
		pszParam = "energy";
		break;

	case eSSminorgroove:				// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		pszParam = "mingrv";
		break;

	case eSSmajorgroove:				// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		pszParam = "majgrv";
		break;

	case eSStwist:					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
		pszParam = "twist";
		break;

	case eSSroll:					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
		pszParam = "roll";
		break;

	case eSStilt:					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
		pszParam = "tilt";
		break;

	case eSSrise:					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
		pszParam = "rise";
		break;

	case eSSslide:					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
		pszParam = "slide";
		break;

	case eSSshift:					// shift int(shift * 10000) e.g 	0.0645 ==> 645
		pszParam = "shift";
		break;

	case eSSrmsd:					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
		pszParam = "rmsd";
		break;

	case eSSORChidVal:
		pszParam = "orchid";
		break;

	default:
		pszParam = "unrecognised";
		break;
	};
return(pszParam);
}




