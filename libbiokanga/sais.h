#pragma once

class CSAIS
{
	void getCounts(const unsigned char *T, int *C, int n, int k, int cs); // find the start or end of each bucket
	void getBuckets(const int *C, int *B, int k, int end);
	void induceSA(const unsigned char *T, int *SA, int *C, int *B, int n, int k, int cs); // compute SA and BWT
	int computeBWT(const unsigned char *T, int *SA, int *C, int *B, int n, int k, int cs);

	/* find the suffix array SA of T[0..n-1] in {0..k-1}^n
		use a working space (excluding T and SA) of at most 2n+O(1) for a constant alphabet */
	int sais_main(const unsigned char *T, int *SA, int fs, int n, int k, int cs, int isbwt);

public:
	CSAIS(void){};
	~CSAIS(void){};

	int sais(const unsigned char *T, int *SA, int n);
	int sais_int(const int *T, int *SA, int n, int k);
	int sais_bwt(const unsigned char *T, unsigned char *U, int *A, int n);
	int sais_int_bwt(const int *T, int *U, int *A, int n, int k);

};
