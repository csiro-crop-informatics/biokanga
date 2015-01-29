// this version has no eBaseEOS processing relying on functionality within uhamming_asm!
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

void
DumpPairA(int SubSeqLen,int Seq1Idx,int Seq2Idx,etSeqBase *pSeq1,etSeqBase *pSeq2)
{
char szSeq1[1000];
char szSeq2[1000];
CSeqTrans::MapSeq2Ascii(pSeq1,SubSeqLen,szSeq1);
CSeqTrans::MapSeq2Ascii(pSeq2,SubSeqLen,szSeq2);
gDiagnostics.DiagOutMsgOnly(eDLInfo,"%1.3d = %s\n\t%1.3d = %s",Seq1Idx,szSeq1,Seq2Idx,szSeq2);
}

// lookup table for quick base complementation
UINT8 MapCpl[] = { eBaseT,		// MapCpl[etBaseA] --> etBaseT
					eBaseG,		// MapCpl[etBaseC] --> etBaseG
					eBaseC,		// MapCpl[etBaseG] --> etBaseC
					eBaseA,		// MapCpl[etBaseT] --> etBaseA
					eBaseN,		// MapCpl[etBaseN] --> etBaseN
					eBaseUndef,	// MapCpl[eBaseUndef] --> eBaseUndef
					eBaseInDel,	// MapCpl[eBaseInDel] --> eBaseInDel
					eBaseEOS};	// MapCpl[eBaseEOS] --> eBaseEOS

// Watson strand only compares only
void GHamDistWatson(UINT8 *pHDs,// where to return Hamming differentials for each subsequence 
			  int SubSeqLen,	// generate Hammings edit distances for subsequences of this length
			  UINT32 SSofs,		// offset between subsequences for current pass
			  UINT8 *pGenomeSeq, // genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
			  UINT32 GenomeLen)	 // total genome length including chrom eBaseEOS markers, but exluding start/end eBaseEOG markers
{
int Idx;

UINT8 PPHD;		// current Hamming for +/+ subsequence edit distances
UINT8 SS1B;		// subsequence 1 5' first base, plus strand
UINT8 SS2B;		// subsequence 2 5' first base, plus strand

UINT8 SE1B;		// subsequence 1 3' last base, plus strand
UINT8 SE2B;		// subsequence 2 3' last base, plus strand

UINT8 *pSS1;	// pts to 5' start of current subsequence1 on plus strand
UINT8 *pSS2;	// pts to 5' start of current subsequence2  on plus strand

UINT8 *pSE1;	// pts to 3' end of current subsequence1 on plus strand
UINT8 *pSE2;	// pts to 3' end of current subsequence2  on plus strand

UINT8 *pHSS1;  // pts to current minimum Huffman distance for subsequence starting at pSS1
UINT8 *pHSS2;  // pts to current minimum Huffman distance for subsequence starting at pSS2

pHSS1 = pHDs;
pHSS2 = &pHDs[SSofs];		

pSS1 = pGenomeSeq;
pSS2 = &pGenomeSeq[SSofs];	

if(SSofs == 41)
	printf("Hello!");

// initial Huffman generation for the subsequences
// this generation is only for the 1st SubSeqLen-1 bases as the final base in subsequences is processed
// within the main iteration loop
// if any eBaseEOS marking chromosome boundaries encountered then the appropriate EOS cnt will be incremented, and
// until these cnts are decremented back to zero in main iteration loop then minimum Hammings will not be updated
PPHD = 0;
pSE1 = pSS1;
pSE2 = pSS2;
for(Idx = 0; Idx < SubSeqLen-1; Idx++)
	{
	if((SE2B = *pSE2++)==eBaseEOG)
		return;
	SE1B = *pSE1++;
	if(SE1B != SE2B)	
		PPHD += 1;
	}

// main iteration loop
// iterate over genome until 3' end of subsequence2 reaches the end of genome marker (eBaseEOS+1)
// end of chromosomes are marked by a eBaseEOS with next chromosome immediately following
// Note: last chrom is terminated by the end of genome marker (eBaseEOS+1) and NOT by eBaseEOS
// 
while((SE2B = *pSE2++) < eBaseEOG)
	{
	SE1B = *pSE1++;
	SS1B = *pSS1++;
	SS2B = *pSS2++;


	// if a difference at 5' of subsequences then Hamming += 1
	// note that Hammings are decremented if difference at 5' subsequences just before end of this iteration
	if(SE1B != SE2B)		
		PPHD += 1;

	if(*pHSS1 > PPHD)		
		*pHSS1 = PPHD;
	if(*pHSS2 > PPHD)		
		*pHSS2 = PPHD;

	pHSS1 += 1;
	pHSS2 += 1;

	// if a difference at 5' start of subsequences then decr Hammings ready for next iteration
	if(SS1B != SS2B)	
		PPHD -= 1;
	}
}



// this version uses the sliding window approach
void 
GHamDistCrick(UINT8 *pHDs,		// where to return Hamming differentials for each subsequence 
			  int SubSeqLen,	// generate Hammings edit distances for subsequences of this length
			  UINT32 SSofs,		// offset between subsequences for current pass
			  UINT8 *pGenomeSeq, // genome sequence (concatenated chrom seqs, separated by eBaseEOSs) with final chrom terminated by eBaseEOG, not eBaseEOS
			  UINT32 GenomeLen)	 // total genome length including chrom eBaseEOS markers, but exluding start/end eBaseEOG markers
{
int Idx;

UINT8 PPHD;		// current Hamming for +/- subsequence edit distances
UINT8 SS1B;		// subsequence 1 5' first base, plus strand
UINT8 SS2B;		// subsequence 2 5' first base, minus strand

UINT8 SE1B;		// subsequence 1 3' last base, plus strand
UINT8 SE2B;		// subsequence 2 3' last base, minus strand

UINT8 *pSS1;	// pts to 5' start of current subsequence1 on plus strand
UINT8 *pSS2;	// pts to 5' start of current subsequence2  on minus strand

UINT8 *pSE1;	// pts to 3' end of current subsequence1 on plus strand
UINT8 *pSE2;	// pts to 3' end of current subsequence2  on minus strand

UINT8 *pHSS1;  // pts to current minimum Huffman distance for subsequence starting at pSS1
UINT8 *pHSS2;  // pts to current minimum Huffman distance for subsequence starting at pSS2

pHSS1 = pHDs;
pHSS2 = &pHDs[GenomeLen - SSofs - SubSeqLen];		

pSS1 = pGenomeSeq;
pSS2 = &pGenomeSeq[GenomeLen - SSofs -1];	


// initial Huffman generation for the subsequences
// this generation is only for the 1st SubSeqLen-1 bases as the final base in subsequences is processed
// within the main iteration loop
PPHD = 0;
pSE1 = pSS1;
pSE2 = pSS2;

for(Idx = 0; Idx < SubSeqLen-1; Idx++)
	{
	if((SE1B = *pSE1++)==eBaseEOG)
		return;
	if((SE2B = *pSE2--)==eBaseEOG)
		return;
	SE2B =MapCpl[SE2B];
	if(SE1B != SE2B)	
		PPHD += 1;
	}

// main iteration loop
// iterate over genome until 3' end of subsequence2 reaches the end of genome marker (eBaseEOS+1)
// end of chromosomes are marked by a eBaseEOS with next chromosome immediately following
// Note: last chrom is terminated by the end of genome marker (eBaseEOS+1) and NOT by eBaseEOS
// 
while((SE2B = *pSE2--) < eBaseEOG)
	{
	SE2B =MapCpl[SE2B];
	SE1B = *pSE1++;
	SS1B = *pSS1++;
	SS2B = MapCpl[*pSS2--];

	// if a difference at 5' of subsequences then Hamming += 1
	// note that Hammings are decremented if difference at 5' subsequences just before end of this iteration
	if(SE1B != SE2B)		
		PPHD += 1;

	if(*pHSS1 > PPHD)		
		*pHSS1 = PPHD;
	if(*pHSS2 > PPHD)		
		*pHSS2 = PPHD;

	pHSS1 += 1;
	pHSS2 -= 1;

	// if a difference at 5' start of subsequences then decr Hammings ready for next iteration
	if(SS1B != SS2B)	
		PPHD -= 1;
	}

if(SSofs == 0)
	return;

// second phase is to wraparound pSS1, and pHSS2, to the end of sequences and continue huffmans after reinitialisation
pHSS2 = &pHDs[GenomeLen - SubSeqLen];		
pSS2 = &pGenomeSeq[GenomeLen - 1];	

// need to re-initial eHuffman generation for the subsequences
// this generation is only for the 1st SubSeqLen-1 bases as the final base in subsequences is processed
// within the main iteration loop
// if any eBaseEOS marking chromosome boundaries encountered then the appropriate EOS cnt will be incremented, and
// until these cnts are decremented back to zero in main iteration loop then minimum Hammings will not be updated
PPHD = 0;
pSE1 = pSS1;
pSE2 = pSS2;
for(Idx = 0; Idx < SubSeqLen-1; Idx++)
	{
	if((SE1B = *pSE1++)==eBaseEOG)
		return;
	if((SE2B = *pSE2--)==eBaseEOG)
		return;
	SE2B = MapCpl[SE2B];
	if(SE1B != SE2B)	
		PPHD += 1;
	}

// main iteration loop
// iterate over genome until 3' end of subsequence2 reaches the end of genome marker (eBaseEOS+1)
// end of chromosomes are marked by a eBaseEOS with next chromosome immediately following
// Note: last chrom is terminated by the end of genome marker (eBaseEOS+1) and NOT by eBaseEOS
// 
while((SE1B = *pSE1++) < eBaseEOG)
	{
	SE2B = MapCpl[*pSE2--];
	SS1B = *pSS1++;
	SS2B = MapCpl[*pSS2--];

	// if a difference at 5' of subsequences then Hamming += 1
	// note that Hammings are decremented if difference at 5' subsequences just before end of this iteration
	if(SE1B != SE2B)		
		PPHD += 1;

	// if subsequences currently do not contain any eBaseEOS then the respective EOS cnt will be zero
	// if EOS cnt is zero then update the minimum Hammings
	if(*pHSS1 > PPHD)		
		*pHSS1 = PPHD;
	if(*pHSS2 > PPHD)
		*pHSS2 = PPHD;

	pHSS1 += 1;
	pHSS2 -= 1;

	// if a difference at 5' start of subsequences then decr Hammings ready for next iteration
	if(SS1B != SS2B)	
		PPHD -= 1;
	}
}
