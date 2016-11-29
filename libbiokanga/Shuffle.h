#pragma once
#include "./commdefs.h"
class CShuffle
{
	CRandomMersenne *m_pRandomMersenne;		// used for generating random numbers larger than RAND_MAX (32767)
	int CHOOSE(int Limit);
	int FChoose(int *p, int N);

public:
	CShuffle(void);
	~CShuffle(void);
	int SeqDPShuffle(int SeqLen,unsigned char *pSeq1, unsigned char *pSeq2);
	int SeqMarkov0(int SeqLen,char *pSeq1, char *pSeq2);
	int SeqMarkov1(int SeqLen, char *pSeq1, char *pSeq2);
	int SeqReverse(int SeqLen, char *pSeq1, char *pSeq2);
	int SeqRegionalShuffle(int SeqLen, char *pSeq1, char *pSeq2, int w);
	int SeqShuffle(int SeqLen, char *pSeq1, char *pSeq2);

};
