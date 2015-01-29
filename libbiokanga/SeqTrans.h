#pragma once
#include "./commdefs.h"

class CSeqTrans
{
public:
	CSeqTrans(void){};
	~CSeqTrans(void){};
	static etSeqBase *MapAscii2Sense(char *pAscii,unsigned int AsciiLen = 0,etSeqBase *pSense=NULL,bool RptMskUpperCase=false);
	static etSeqBase *MapAscii2Antisense(char *pAscii,unsigned int AsciiLen = 0,etSeqBase *pAntisense=NULL,bool RptMskUpperCase=false);
	static void ReverseComplement(unsigned int SeqLen,etSeqBase *pSeq);
	static char * MapSeq2Ascii(etSeqBase *pBases,unsigned int SeqLen,char *pAscii = NULL,
							char BaseN = 'N',				// what to map eBaseN onto
							char BaseUndef = 'U',			// what to map eBaseUndef onto
							char BaseInDel = '-',			// what to map eBaseInDel onto
							bool RptMskUpperCase=false);

	static char *MapPackedSeq2Ascii(int StartNibble,// which nibble to start from - 0 if lower or first nibble
							etSeqBase *pBases,		// packed bases (bases in lower and upper nibbles) to map to ascii
							unsigned int SeqLen,	// total number of bases to map
							char *pAscii = NULL,	// where to return ascii - NOTE: will have terminating '\0' appended
							char BaseN = 'N',				// what to map eBaseN onto
							char BaseUndef = 'U',			// what to map eBaseUndef onto
							char BaseInDel = '-',			// what to map eBaseInDel onto
							bool RptMskUpperCase=false);  // if true then uppercase represents softmasked repeats

	static char MapBase2Ascii(etSeqBase Base,				// base to map 
							char BaseN = 'N',				// what to map eBaseN onto
							char BaseUndef = 'U',			// what to map eBaseUndef onto
							char BaseInDel = '-',			// what to map eBaseInDel onto
							bool RptMskUpperCase=false);

	static char *MapSeq2UCAscii(etSeqBase *pBases,			// bases to map to uppercase ascii
							unsigned int SeqLen,			// number of bases to map
							char *pAscii,					// where to return ascii - NOTE: will have terminating '\0' appended
							char BaseN = 'N',				// what to map eBaseN onto
							char BaseUndef = 'U',				// what to map eBaseUndef onto
							char BaseInDel = '-');				// what to map eBaseInDel onto

	static char *MapSeq2LCAscii(etSeqBase *pBases,			// bases to map to lowercase ascii
							unsigned int SeqLen,			// number of bases to map
							char *pAscii,					// where to return ascii - NOTE: will have terminating '\0' appended
							char BaseN = 'n',				// what to map eBaseN onto
							char BaseUndef = 'u',				// what to map eBaseUndef onto
							char BaseInDel = '-');				// what to map eBaseInDel onto

	static void ReverseSeq(unsigned int SeqLen,etSeqBase *pSeq);
	static	void ComplementStrand(unsigned int SeqLen,etSeqBase *pSeq);
	static void RemoveMasking(etSeqBase *pBases,int SeqLen);
};
