/*
 * sais.c for sais-lite
 * Copyright (c) 2008-2009 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
// Changes made to original source code are purely to make function scope within the CSAIS class framework
#include "stdafx.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if _WIN32
#include <process.h>
#include "./commhdrs.h"
#else
#include <sys/mman.h>
#include <pthread.h>
#include "./commhdrs.h"
#endif

//#include <stdlib.h>
#ifdef _OPENMP
# include <omp.h>
#endif
#include "sais.h"


#define chr(i) (cs == sizeof(int) ? ((const int *)T)[i]:((const unsigned char *)T)[i])

/* find the start or end of each bucket */
void
CSAIS::getCounts(const unsigned char *T, int *C, int n, int k, int cs) {
#ifdef _OPENMP
  int *D;
  int i, j, p, sum, first, last;
  int thnum, maxthreads = omp_get_max_threads();
#pragma omp parallel default(shared) private(D, i, thnum, first, last)
  {
    thnum = omp_get_thread_num();
    D = C + thnum * k;
    first = n / maxthreads * thnum;
    last = (thnum < (maxthreads - 1)) ? n / maxthreads * (thnum + 1) : n;
    for(i = 0; i < k; ++i) { D[i] = 0; }
    for(i = first; i < last; ++i) { ++D[chr(i)]; }
  }
  if(1 < maxthreads) {
#pragma omp parallel for default(shared) private(i, j, p, sum)
    for(i = 0; i < k; ++i) {
      for(j = 1, p = i + k, sum = C[i]; j < maxthreads; ++j, p += k) {
        sum += C[p];
      }
      C[i] = sum;
    }
  }
#else
  int i;
  for(i = 0; i < k; ++i) { C[i] = 0; }
  for(i = 0; i < n; ++i) { ++C[chr(i)]; }
#endif
}

void
CSAIS::getBuckets(const int *C, int *B, int k, int end) {
  int i, sum = 0;
  if(end) { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum; } }
  else { for(i = 0; i < k; ++i) { sum += C[i]; B[i] = sum - C[i]; } }
}

/* compute SA and BWT */
void
CSAIS::induceSA(const unsigned char *T, int *SA, int *C, int *B, int n, int k, int cs) {
  int *b, i, j;
  int c0, c1;
  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
    j = SA[i], SA[i] = ~j;
    if(0 < j) {
      --j;
      if((c0 = chr(j)) != c1) 
		{ 
		B[c1] = (int)(b - SA); 
		b = SA + B[c1 = c0]; 
	    }
      *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
    }
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      --j;
      if((c0 = chr(j)) != c1) 
		{ 
		B[c1] = (int)(b - SA); 
		b = SA + B[c1 = c0]; 
		}
      *--b = ((j == 0) || (chr(j - 1) > c1)) ? ~j : j;
    } else {
      SA[i] = ~j;
    }
  }
}

int
CSAIS::computeBWT(const unsigned char *T, int *SA, int *C, int *B, int n, int k, int cs) {
  int *b, i, j, pidx = -1;
  int c0, c1;
  /* compute SAl */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 0); /* find starts of buckets */
  j = n - 1;
  b = SA + B[c1 = chr(j)];
  *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
  for(i = 0; i < n; ++i) {
    if(0 < (j = SA[i])) {
      --j;
      SA[i] = ~(c0 = chr(j));
      if(c0 != c1) 
		{ 
		B[c1] = (int)(b - SA); 
		b = SA + B[c1 = c0]; 
		}
      *b++ = ((0 < j) && (chr(j - 1) < c1)) ? ~j : j;
    } else if(j != 0) {
      SA[i] = ~j;
    }
  }
  /* compute SAs */
  if(C == B) { getCounts(T, C, n, k, cs); }
  getBuckets(C, B, k, 1); /* find ends of buckets */
  for(i = n - 1, b = SA + B[c1 = 0]; 0 <= i; --i) {
    if(0 < (j = SA[i])) {
      --j;
      SA[i] = (c0 = chr(j));
      if(c0 != c1) 
		{ 
		B[c1] = (int)(b - SA); 
		b = SA + B[c1 = c0]; 
		}
      *--b = ((0 < j) && (chr(j - 1) > c1)) ? ~((int)chr(j - 1)) : j;
    } else if(j != 0) {
      SA[i] = ~j;
    } else {
      pidx = i;
    }
  }
  return pidx;
}


/* find the suffix array SA of T[0..n-1] in {0..k-1}^n
   use a working space (excluding T and SA) of at most 2n+O(1) for a constant alphabet */
int
CSAIS::sais_main(const unsigned char *T, int *SA, int fs, int n, int k, int cs, int isbwt) {
  int *C, *B, *RA;
  int i, j, c, m, p, q, plen, qlen, name, pidx = 0;
  int c0, c1;
  int diff;
#ifdef _OPENMP
  int maxthreads = omp_get_max_threads();
#else
# define maxthreads 1
#endif

  /* stage 1: reduce the problem by at least 1/2
     sort all the S-substrings */
  if((maxthreads * k) <= fs) {
    C = SA + n;
    B = ((1 < maxthreads) || (k <= (fs - k))) ? C + k : C;
  } else {
    if((C = (int *)malloc(maxthreads * k * sizeof(int))) == NULL) { return -2; }
    B = (1 < maxthreads) ? C + k : C;
  }
  getCounts(T, C, n, k, cs); getBuckets(C, B, k, 1); /* find ends of buckets */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
  for(i = 0; i < n; ++i) { SA[i] = 0; }
  for(i = n - 2, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
    if((c0 = chr(i)) < (c1 + c)) { c = 1; }
    else if(c != 0) { SA[--B[c1]] = i + 1, c = 0; }
  }
  induceSA(T, SA, C, B, n, k, cs);
  if(fs < (maxthreads * k)) { free(C); }

  /* compact all the sorted substrings into the first m items of SA
     2*m must be not larger than n (proveable) */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i, j, p, c0, c1)
  for(i = 0; i < n; ++i) {
    p = SA[i];
    if((0 < p) && (chr(p - 1) > (c0 = chr(p)))) {
      for(j = p + 1; (j < n) && (c0 == (c1 = chr(j))); ++j) { }
      if((j < n) && (c0 < c1)) { SA[i] = ~p; }
    }
  }
  for(i = 0, m = 0; i < n; ++i) { if((p = SA[i]) < 0) { SA[m++] = ~p; } }
#else
  for(i = 0, m = 0; i < n; ++i) {
    p = SA[i];
    if((0 < p) && (chr(p - 1) > (c0 = chr(p)))) {
      for(j = p + 1; (j < n) && (c0 == (c1 = chr(j))); ++j) { }
      if((j < n) && (c0 < c1)) { SA[m++] = p; }
    }
  }
#endif
  j = m + (n >> 1);
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
  for(i = m; i < j; ++i) { SA[i] = 0; } /* init the name array buffer */
  /* store the length of all substrings */
  for(i = n - 2, j = n, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
    if((c0 = chr(i)) < (c1 + c)) { c = 1; }
    else if(c != 0) { SA[m + ((i + 1) >> 1)] = j - i - 1; j = i + 1; c = 0; }
  }
  /* find the lexicographic names of all substrings */
  for(i = 0, name = 0, q = n, qlen = 0; i < m; ++i) {
    p = SA[i], plen = SA[m + (p >> 1)], diff = 1;
    if(plen == qlen) {
      for(j = 0; (j < plen) && (chr(p + j) == chr(q + j)); ++j) { }
      if(j == plen) { diff = 0; }
    }
    if(diff != 0) { ++name, q = p, qlen = plen; }
    SA[m + (p >> 1)] = name;
  }

  /* stage 2: solve the reduced problem
     recurse if names are not yet unique */
  if(name < m) {
    RA = SA + n + fs - m;
    for(i = m + (n >> 1) - 1, j = m - 1; m <= i; --i) {
      if(SA[i] != 0) { RA[j--] = SA[i] - 1; }
    }
    if(sais_main((unsigned char *)RA, SA, fs + n - m * 2, m, name, sizeof(int), 0) != 0) { return -2; }
    for(i = n - 2, j = m - 1, c = 0, c1 = chr(n - 1); 0 <= i; --i, c1 = c0) {
      if((c0 = chr(i)) < (c1 + c)) { c = 1; }
      else if(c != 0) { RA[j--] = i + 1, c = 0; } /* get p1 */
    }
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for(i = 0; i < m; ++i) { SA[i] = RA[SA[i]]; } /* get index */
  }

  /* stage 3: induce the result for the original problem */
  if((maxthreads * k) <= fs) {
    C = SA + n;
    B = ((1 < maxthreads) || (k <= (fs - k))) ? C + k : C;
  } else {
    if((C = (int *)malloc(maxthreads * k * sizeof(int))) == NULL) { return -2; }
    B = (1 < maxthreads) ? C + k : C;
  }
  /* put all left-most S characters into their buckets */
  getCounts(T, C, n, k, cs); getBuckets(C, B, k, 1); /* find ends of buckets */
#ifdef _OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
  for(i = m; i < n; ++i) { SA[i] = 0; } /* init SA[m..n-1] */
  for(i = m - 1; 0 <= i; --i) {
    j = SA[i], SA[i] = 0;
    SA[--B[chr(j)]] = j;
  }
  if(isbwt == 0) { induceSA(T, SA, C, B, n, k, cs); }
  else { pidx = computeBWT(T, SA, C, B, n, k, cs); }
  if(fs < (maxthreads * k)) { free(C); }

  return pidx;
#ifndef _OPENMP
# undef maxthreads
#endif
}

int
CSAIS::sais(const unsigned char *T, int *SA, int n) {
  if((T == NULL) || (SA == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { SA[0] = 0; } return 0; }
  return sais_main(T, SA, 0, n, 256, sizeof(unsigned char), 0);
}

int
CSAIS::sais_int(const int *T, int *SA, int n, int k) {
  if((T == NULL) || (SA == NULL) || (n < 0) || (k <= 0)) { return -1; }
  if(n <= 1) { if(n == 1) { SA[0] = 0; } return 0; }
  return sais_main((const unsigned char *)T, SA, 0, n, k, sizeof(int), 0);
}

int
CSAIS::sais_bwt(const unsigned char *T, unsigned char *U, int *A, int n) {
  int i, pidx;
  if((T == NULL) || (U == NULL) || (A == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }
  pidx = sais_main(T, A, 0, n, 256, sizeof(unsigned char), 1);
  if(pidx < 0) { return pidx; }
  U[0] = T[n - 1];
  for(i = 0; i < pidx; ++i) { U[i + 1] = (unsigned char)A[i]; }
  for(i += 1; i < n; ++i) { U[i] = (unsigned char)A[i]; }
  pidx += 1;
  return pidx;
}

int
CSAIS::sais_int_bwt(const int *T, int *U, int *A, int n, int k) {
  int i, pidx;
  if((T == NULL) || (U == NULL) || (A == NULL) || (n < 0)) { return -1; }
  if(n <= 1) { if(n == 1) { U[0] = T[0]; } return n; }
  pidx = sais_main((const unsigned char *)T, A, 0, n, k, sizeof(int), 1);
  if(pidx < 0) { return pidx; }
  U[0] = T[n - 1];
  for(i = 0; i < pidx; ++i) { U[i + 1] = A[i]; }
  for(i += 1; i < n; ++i) { U[i] = A[i]; }
  pidx += 1;
  return pidx;
}
