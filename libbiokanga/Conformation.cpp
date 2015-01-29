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

CConformation::CConformation(void)
{
m_pStructParams = NULL;
}

CConformation::~CConformation(void)
{
if(m_pStructParams != NULL)
	delete m_pStructParams;
m_pStructParams = NULL;
}

// StructParamsLoaded
// Returns true if structural parameters loaded
bool
CConformation::StructParamsLoaded(void)
{
return(m_pStructParams == NULL ? false : true);
}

// LoadStructParams
// Load structural parameters from file
// File layout is that of one comma separated set of parameters per line with optional heading line
// Octamer,Twist,Roll,Tilt,Rise,Slide,Shift,3-StepTwist,3-StepRoll,3-StepSlide,3-StepShift,Energy,
// MinorGroove,RMSD,Q-Twist,Q+Twist,Q-Roll,Q+Roll,3Q-Twist,3Q+Twist,3Q-Roll,3Q+Roll, ORChID (hydroxyl radical cleavage values)
//
teBSFrsltCodes 
CConformation::LoadStructParams(char *pszStructParamsFile) // load structural parameters, from internal fixed constant array or from file
{
FILE *pParamsStream;
int LineNum;
int NumOctParams;
char szLineBuff[512];
int Cnt;
char Octamer[9];
etSeqBase octbases[8];
int OctIdx;
double Twist,Roll,Tilt,Rise,Slide,Shift,TriStepTwist,TriStepRoll,TriStepSlide,TriStepShift;
double Energy,MinorGroove,RMSD,QminusTwist,QplusTwist,QminusRoll,QplusRoll;
double TriQminusTwist,TriQplusTwist,TriQminusRoll,TriQplusRoll,ORChID;

tsStructParam *pStruct1;
char chr, *pDst, *pSrc;

if(m_pStructParams == NULL)
	{	
	m_pStructParams = (tsStructParam *) new unsigned char[cStructParamAllocSize];
	if(m_pStructParams == NULL)
		{
		AddErrMsg("CConformation::LoadStructParams","Unable to allocate memory to hold structural parameters");
		return(eBSFerrMem);
		}
	memset(m_pStructParams,0,cStructParamAllocSize);
	}
m_szStructParamFile[0] = '\0';


if((pParamsStream = fopen(pszStructParamsFile,"r"))==NULL)
	{
	AddErrMsg("CConformation::LoadStructParams","Unable to open parameters file %s error: %s",pszStructParamsFile,strerror(errno));
	delete m_pStructParams;
	m_pStructParams = NULL;
	return(eBSFerrOpnFile);
	}


LineNum = 0;
NumOctParams = 0;
int NxScanPsn;
int ParamLines = 0;
while(fgets(szLineBuff,sizeof(szLineBuff),pParamsStream)!= NULL)
	{
	LineNum++;
	if(strlen(szLineBuff) < 5)	// simply slough lines which are too short to contain anything worth parsing
		continue;

	// strip any whitespace and quotes
	pDst = pSrc = szLineBuff;
	while(chr = *pSrc++)
		if(!isspace(chr) && chr != '\'' && chr != '"')
			*pDst++ = chr;
	*pDst = '\0';
	if(szLineBuff[0] == '\0')
		continue;
	ParamLines += 1;
	// if unable to initially parse the first line then assume it's a header line and slough
	 Cnt = sscanf(szLineBuff,"%8[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf %n",
			Octamer,&Twist,&Roll,&Tilt,&Rise,&Slide,&Shift,&TriStepTwist,&TriStepRoll,&TriStepSlide,&TriStepShift, &NxScanPsn);
	 if(Cnt < 11 && ParamLines == 1)
		 continue;

	 if(Cnt == 11)
		 Cnt += sscanf(&szLineBuff[NxScanPsn],",%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
				&Energy,&MinorGroove,&RMSD,&QminusTwist,&QplusTwist,&QminusRoll,&QplusRoll,&TriQminusTwist,&TriQplusTwist,&TriQminusRoll,&TriQplusRoll, &ORChID);

	 if(Cnt != 23)
		{
		AddErrMsg("CConformation::LoadStructParams","Error parsing structural parameters file %s at line %d, expected 23 but only parsed %d parameters\n%s\n",pszStructParamsFile,LineNum,Cnt,szLineBuff);
		fclose(pParamsStream);
		delete m_pStructParams;
		m_pStructParams = NULL;
		return(eBSFerrStructParm); 
		}

	CSeqTrans::MapAscii2Sense(Octamer,8,octbases);
	if((OctIdx = StructParamIdx(octbases)) < 0)
		{
		AddErrMsg("CConformation::LoadStructParams","Octamer in file %s at line %d contains unrecognised bases (only 'a','c','g','t' accepted)\n%s\n",pszStructParamsFile,LineNum,szLineBuff);
		fclose(pParamsStream);
		delete m_pStructParams;
		m_pStructParams = NULL;
		return((teBSFrsltCodes)OctIdx); 
		}
	pStruct1	 = &m_pStructParams[OctIdx];
	if(pStruct1->Param.twist != 0)
		{
		fclose(pParamsStream);
		delete m_pStructParams;
		AddErrMsg("CConformation::LoadStructParams","Duplicate octamer structural properties parsed at line %d in file %s\n",LineNum,pszStructParamsFile);
		return(eBSFerrStructParm);
		}

	pStruct1->Param.energy = (int)(Energy * 10000.0);
	pStruct1->Param.minorgroove = (int)(MinorGroove * 10000.0);
	pStruct1->Param.twist =  (int)(Twist * 10000.0);
	pStruct1->Param.roll = (int)(Roll * 10000.0);
	pStruct1->Param.tilt = (int)(Tilt * 10000.0);
	pStruct1->Param.rise = (int)(Rise * 10000.0);
	pStruct1->Param.slide = (int)(Slide * 10000.0);
	pStruct1->Param.shift = (int)(Shift * 10000.0);
	pStruct1->Param.rmsd = (int)(RMSD * 10000.0);
	pStruct1->Param.orchid = (int)(ORChID * 10000.0);
	NumOctParams++;
	}
fclose(pParamsStream);

if(NumOctParams != cNumParamOctamers)
	AddErrMsg("CConformation::LoadStructParams","Warning, not all octamers in '%s' have structural properties\n",pszStructParamsFile);
strncpy(m_szStructParamFile,pszStructParamsFile,_MAX_PATH-1);
m_szStructParamFile[_MAX_PATH-1] = '\0';
return(eBSFSuccess);
}

// reproducably psudeo-randomise conformation characteristics
//
int
CConformation::PsudeoRandomise(void)
{
int Len;
int OctIdx;
tsStructParam tmp;
TRandomMersenne Random(1); // random sequence must be reproducable so seed required
for(Len = cNumParamOctamers; Len > 1; Len--)
	{				
    OctIdx       = Random.IRandom(0,Len-1);
    tmp         = m_pStructParams[OctIdx];
    m_pStructParams[OctIdx]   = m_pStructParams[Len-1];
    m_pStructParams[Len-1] = tmp;
	}
return(eBSFSuccess);
}


// determine idex to use for assumed octamer sequence
// eBSFerrStructStep returned if any base in octamer is indeterminate - 'N'
int
CConformation::StructParamIdx(etSeqBase *pOctamer)		// sequence
{
etSeqBase Base;
int OctIdx = 0;
int Len = 8;
if(pOctamer == NULL || m_pStructParams == NULL)
	return((int)eBSFerrParams);

while(Len--)
	{
	OctIdx <<= 2;
	Base = *pOctamer++  & ~cRptMskFlg;
	if(Base > eBaseT)
		return((int)eBSFerrStructStep);								// unrecognised base
	OctIdx |= Base;
	}
return(OctIdx);
}

teBSFrsltCodes
CConformation::GetSequenceConformation(teStructStats Param,		// which structural parameter value to return
				 unsigned int iStartOfs, // initial starting offset (0..n) in pSeq
				  unsigned int iNumSteps,		  // number of steps (0 for all) to process starting at pSeq[iStartPsn]|pSeq[iStartPsn+1]
  				  unsigned int SeqLen,			  // total length of sequence
				  etSeqBase *pSeq,				  // sequence to be processed
				  int *pRetConfValue,			  // where to return conformation
				  int UndefBaseValue)			  // value to return for undefined or indeterminate ('N') bases 
{
unsigned int Step;
unsigned int LastStep;
int IdxRetVal;

if(SeqLen < 8 || iStartOfs >= SeqLen - 1 || 
   iStartOfs + iNumSteps >= SeqLen ||
   pSeq == NULL || m_pStructParams == NULL)
	return(eBSFerrParams);

if(iNumSteps == 0)
	iNumSteps = SeqLen - iStartOfs - 1;
LastStep = iStartOfs + iNumSteps;

for(IdxRetVal =0,Step = iStartOfs+1; Step <= LastStep; IdxRetVal++,Step++)
	pRetConfValue[IdxRetVal] = StructValue(Param,Step,SeqLen,pSeq,UndefBaseValue);
return(eBSFSuccess);
}

int
CConformation::StructValue(teStructStats Param,		// which structural parameter value to return
			unsigned int Step,			// which step in sequence to return structural value for
			unsigned int SeqLen,		// total length of sequence
			etSeqBase *pSeq,			// sequence to be processed
			int UndefBaseValue)			// value to return for undefined or indeterminate ('N') bases 
{
tsStructParam *pStruct;
etSeqBase *pOctamer;
etSeqBase octamer[8];
int OctIdx;
unsigned int InterpStep;
int ParamOfs;
bool bContinue;
bool bLeft;
int Iters;
INT64 InterpValue;
bContinue = false;
int	Cnt;
int	OctOfs;
int	SeqOfs;

switch(Param) {
	case eSSenergy:				// minimal energy int(energy * 10000) e.g. -408.2632 ==> -4082632
		ParamOfs = offsetof(tsStructParam,Param.energy);
		break;

	case eSSminorgroove:				// groove int(dimensions * 10000) e.g 10.784 ==> 107840
		ParamOfs = offsetof(tsStructParam,Param.minorgroove);
		break;

	case eSStwist:					// twist int(angle * 10000) e.g 	37.6262 ==> 376262
		ParamOfs = offsetof(tsStructParam,Param.twist);
		break;

	case eSSroll:					// roll int(angle * 10000) e.g 	2.4139 ==> 24139
		ParamOfs = offsetof(tsStructParam,Param.roll);
		break;

	case eSStilt:					// tilt int(tilt * 10000) e.g 	-0.022 ==> -220
		ParamOfs = offsetof(tsStructParam,Param.tilt);
		break;

	case eSSrise:					// rise int(rise * 10000) e.g 	3.1409 ==> 31409
		ParamOfs = offsetof(tsStructParam,Param.rise);
		break;

	case eSSslide:					// slide int(slide * 10000) e.g 	-0.0968 ==> -968	
		ParamOfs = offsetof(tsStructParam,Param.slide);
		break;

	case eSSshift:					// shift int(shift * 10000) e.g 	0.0645 ==> 645
		ParamOfs = offsetof(tsStructParam,Param.shift);
		break;

	case eSSrmsd:					// rmsd int(rmsd * 10000) e.g 	0.3078 ==> 3078
		ParamOfs = offsetof(tsStructParam,Param.rmsd);
		break;

	case eSSORChidVal:
		ParamOfs = offsetof(tsStructParam,Param.orchid);
		break;

	default:
		return(eBSFerrParams);
	};

// if can use the midstep of an octamer then do so...
if(Step > 3 && Step < SeqLen - 3)
	{
	pOctamer = pSeq + Step - 4;
	if((OctIdx = StructParamIdx(pOctamer)) < 0)
		return(UndefBaseValue);
	else
		{
		pStruct = &m_pStructParams[OctIdx];
		return(*(int *)(((unsigned char *)pStruct)+ParamOfs));
		}
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
InterpStep = Step;
memset(octamer,eBaseA,8);
if(InterpStep < 4)
	{
	Cnt = 4 + InterpStep;
	OctOfs = 4 - InterpStep;
	SeqOfs = 0;
	}
else
	{
	InterpStep = 8 - (SeqLen - InterpStep);						// adjust Step to 5..7
	Cnt = 12 - InterpStep;
	OctOfs = 0;
	SeqOfs = SeqLen - (Cnt + 1);
	}
memmove(&octamer[OctOfs],&pSeq[SeqOfs],Cnt);

//		Interpolate(Step,octamer,pRetStructParam);
if(InterpStep > 4)								// sequence to left is known
	{
	bLeft = true;
	InterpStep -= 4;						
	}
else										// else sequence to right is known
	{
	bLeft = false;
	InterpStep = 4 - InterpStep;	
	}

switch(InterpStep) {
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

InterpValue=0;
for(Cnt = 0; Cnt < Iters; Cnt++)
	{
	OctIdx = StructParamIdx(octamer);
	if(OctIdx < 0)
		return(UndefBaseValue);
	pStruct = &m_pStructParams[OctIdx];
	InterpValue += *(int *)(((unsigned char *)pStruct)+ParamOfs);;

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
InterpValue /= Iters;
return((int)InterpValue);
}



	