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

CStats::CStats(void)
{
m_pLogFact = NULL;		// holds array of all factorials from 1..m_AllocLogFacts as log(fact)
Init();
}

CStats::~CStats(void)
{
if(m_pLogFact != NULL)
	delete m_pLogFact;
}

void
CStats::Init(void)
{
if(m_pLogFact != NULL)
	{
	delete m_pLogFact;
	m_pLogFact = NULL;
	}
m_NumLogFacts = 0;		// largest precalculated log(factorial) in m_pLogFact
m_AllocLogFacts = 0;	// how many entries have been allocated in m_pLogFact
m_pLogFact = NULL;		// holds array of all factorials from 1..m_AllocLogFacts as log(fact)
m_lncof[0] = 76.18009173;
m_lncof[1] = -86.50532033;
m_lncof[2] = 24.01409822; 
m_lncof[3] = -1.231739516; 
m_lncof[4] = 0.00120858003;
m_lncof[5] = -0.536382e-5;
}

/* Returns log(Gamma(z))  for  */
/*   Gamma(a) = int_0^infty x^{a-1} exp(-x) dx */
#define DTOL 1e-8  /* For convergence */
double 
CStats::lngam(double z)
{ 
double retval;

if(z<=0.0)
	return(-1.0);

if (z<1.0)  
	retval=lngam(z+1)-log(z);
else
	{ 
	int i;
	double zz, ser, ztmp = z + 4.50 - (z-0.5)*log(z + 4.50);
	i=0;  zz=z-1;  ser=1.0;
	while(i<=5)  
		ser += m_lncof[i++]/(zz += 1.0);
	retval = -ztmp + log(2.50662827465*ser); 
	}
return retval; 
}

/* Workhorse function: */
/* Incomplete gamma function with scale=1: */
/* Returns  Prob(Gam(alpha,1)<=xx)  */
/*   = (1/Ga(alpha)) int(0,xx) y^{alpha-1} exp(-y) dy  */
/* Adapted from Press etal, 2nd ed, p216-219. */

double 
CStats::gamminc(double aa, double xx)
{ 
if (aa<=0.0)
	return(-1);
if (xx<=0.0)  
	return 0.0;
if (xx < aa+1) /* Use power-series expansion */
	{ 
	double ap=aa, sum, term;  term=sum=1.0/ap;
	do  
		sum += (term *= xx/(ap += 1.0));
	while (fabs(term) >= fabs(sum)*DTOL);
	return(sum*exp(-xx + aa*log(xx) - lngam(aa)));
	}
else /* Use continued-fraction expansion */
    { 
	int n;  double gold=0.0, fac=1.0;
    double a0=1.0, a1=xx, b0=0.0, b1=1.0;
    for (n=1; ; n++)
		{ 
		double anf=n*fac, ana=n-aa;
        a0=(a1+a0*ana)*fac;   b0=(b1+b0*ana)*fac;
        a1=xx*a0+anf*a1;      b1=xx*b0+anf*b1;
        if (fabs(a1) > 1e-11)
			{
			double g=b1*(fac=1.0/a1);
            if (fabs(g-gold) < DTOL)
              return  1.0 - g*exp(-xx + aa*log(xx) - lngam(aa));
            gold=g; 
			}
		}
	}
}



/* Chi-square P-value  */
/* Returns  Prob( Chi(df)>=xx )   */
/* df is the degrees of freedom for a chi-square distribution */

double 
CStats::ChiSqr2PVal(int df, double ChiSqr)
{ 
if (df < 1 || ChiSqr<0.0)
   return(-1);
return (1.0 - gamminc(((double)df)/2, ChiSqr/2));
}


// Cells expected to be in row order, row 1 containing Cols values is followed by row 2... row Rows.
// Thus for a 2*2 contingency table:
// Rows == 2
// Cols == 2
// pCell[0] == sample1 true count
// pCell[1] == sample1 false count
// pCell[2] == sample2 true count
// pCell[3] == sample2 false count
//
// Returns Ch-square with Yates Continuity correction
// If any expected value is less than 5 then returns -1 as error
double
CStats::CalcChiSqr(int Rows,int Cols,int *pCells)
{
int IdxJ,IdxI;
int *pCell;
double *pRowTotal;
double *pColTotal;
double *pExpect;
double RowItotals[cMaxChiSqrRows];
double ColJtotals[cMaxChiSqrCols];
double Expects[cMaxChiSqrRows][cMaxChiSqrCols];
double TotalF;
double ChiSqr;
double tmp;
memset(RowItotals,0,sizeof(double)*Rows);
memset(ColJtotals,0,sizeof(double)*Cols);

// total rows and columns-
pCell = pCells;
pRowTotal = RowItotals;
TotalF = 0.0;
for(IdxI = 0; IdxI < Rows; IdxI++,pRowTotal++)
	{
	pColTotal = ColJtotals;
	for(IdxJ=0;IdxJ<Cols;IdxJ++,pCell++,pColTotal++)
		{
		TotalF += *pCell;
		*pRowTotal += *pCell;
		*pColTotal += *pCell;
		}
	}
// totals now known, calc expectation for each cell
pExpect = Expects[0];
pRowTotal = RowItotals;
for(IdxI = 0; IdxI < Rows; IdxI++,pRowTotal++)
	{
	pColTotal = ColJtotals;
	for(IdxJ=0;IdxJ<Cols;IdxJ++,pExpect++,pColTotal++)
		{
		*pExpect = (*pColTotal * (*pRowTotal)) / TotalF;
		if(*pExpect < 5.0)		
			return(-1);
		}
	}

// Expects now known, calc ChiSqr
pExpect = Expects[0];
pCell = pCells;
pRowTotal = RowItotals;
ChiSqr = 0.0;
for(IdxI = 0; IdxI < Rows; IdxI++,pRowTotal++)
	{
	pColTotal = ColJtotals;
	for(IdxJ=0;IdxJ<Cols;IdxJ++,pExpect++,pCell++,pColTotal++)
		{
		if(*pExpect > 0.0)
			{
			tmp = fabs(((double)*pCell) - *pExpect)-0.5; // yates continuity correction
			ChiSqr += (tmp * tmp) / *pExpect;
			}
		}
	}
return(ChiSqr);
}


// FishersExactTest
// Returns one tailed P value
// P0 = (R1T! * R2T! * C1T! * C2T!) / (TotalN! * R1C1! * R1C2! * R2C1! * R2C2!)
// If returned value is negative then: 
// -1.0 at least one cell value was < 0
// -2.0 unable to allocate required memory
//
double
CStats::FishersExactTest(int R1C1,		// sample1 true
						 int R1C2,		// sample1 false
						 int R2C1,		// sample2 true
						 int R2C2)		// sample2 false
{
double Num,PValue;
int R1T,R2T,C1T,C2T,TotalN;
double *pLogFact;
INT64 Row1Tot;
INT64 Row2Tot;
INT64 TargCnt;int tmp1,tmp2;
int Idx;

// sanity check
if(R1C1 < 0 || R1C2 < 0 || R2C1 < 0 || R2C2 < 0)
	return(-1.0);

// scale samples so that total counts are less than cMaxTotSampleCnt, scaling is per sample row
// otherwise likely to hit memory resource limits for holding precalc'd log factorials
// any way, if user has sample counts of this magnitude then why not use chi-square?

while((R1C1 + R1C2 + R2C1 + R2C2) > cMaxTotSampleCnt)
	{
	Row1Tot = R1C1 + R1C2;
	Row2Tot = R2C1 + R2C2;

	if(Row1Tot > Row2Tot)
		{
		if(Row2Tot < (cMaxTotSampleCnt / 2))
			TargCnt = cMaxTotSampleCnt - Row2Tot;
		else
			TargCnt = cMaxTotSampleCnt / 2;
		R1C1 = (int)(((INT64)R1C1 * TargCnt)/Row1Tot);
		R1C2 = (int)(((INT64)R1C2 * TargCnt)/Row1Tot);
		}
	else
		{
		if(Row1Tot < (cMaxTotSampleCnt / 2))
			TargCnt = cMaxTotSampleCnt - Row1Tot;
		else
			TargCnt = cMaxTotSampleCnt / 2;
		R2C1 = (int)(((INT64)R2C1 * TargCnt)/Row2Tot);
		R2C2 = (int)(((INT64)R2C2 * TargCnt)/Row2Tot);
		}
	}

// reorder so that lowest valued is in R1C1 by swaping rows and columns
// this is worth the cost as it minimises that number of iterations required when summing P values later
if((R1C1 > R1C2 &&  R1C2 < R2C1) || (R1C1 > R2C2 && R2C2 < R2C1))	// need to exchange cols?
	{
	tmp1 = R1C1;
	tmp2 = R2C1;
	R1C1 = R1C2;
	R2C1 = R2C2;
	R1C2 = tmp1;
	R2C2 = tmp2;
	}

if(R1C1 > R2C1)						// need to exchange rows?
	{
	tmp1 = R1C1;
	tmp2 = R1C2;
	R1C1 = R2C1;
	R1C2 = R2C2;
	R2C1 = tmp1;
	R2C2 = tmp2;
	}
// now ordered such that R1C1 is lowest valued

// to improve performance then scale counts down such that the minimum value is clamped to be no more than 1000000
// surely with minimum counts of 1000000+ then user should be using chi-square
if(R1C1 > 100000)
	{
	R1C2 = (int)(((INT64)1000000 * R1C2)/R1C1);
	R2C1 = (int)(((INT64)1000000 * R2C1)/R1C1);
	R2C2 = (int)(((INT64)1000000 * R2C2)/R1C1);
	R1C1 = 1000000;
	}

// determine row,col, plus total sums
R1T = R1C1 + R1C2;
R2T = R2C1 + R2C2;
C1T = R1C1 + R2C1;
C2T = R1C2 + R2C2;
TotalN = R1T + R2T;

// check if the precalc'd log factorials can handle TotalN!
// if not then need to realloc to hold additional factorials
if(m_pLogFact == NULL || TotalN >= m_NumLogFacts)
	{
	try {
		pLogFact = (double *)new double [TotalN + cAllocLogFacts];
		}
	catch(...) {
		// assume alloc failed because already using too much memory
		delete m_pLogFact;
		m_pLogFact = NULL;
		m_AllocLogFacts = 0;
		m_NumLogFacts = 0;
		pLogFact = NULL;
		}
	if(pLogFact == NULL)	// will be NULL if prev failed, released memory so try again
		pLogFact = (double *)new double [TotalN + cAllocLogFacts];
	
	if(pLogFact == NULL)
		return(-2.0);

	if(m_pLogFact != NULL)
		{
		if(m_NumLogFacts)
			memmove(pLogFact,m_pLogFact,sizeof(double) * m_NumLogFacts);
		delete m_pLogFact;
		}
	else
		m_NumLogFacts = 0;
	m_pLogFact = pLogFact;
	m_AllocLogFacts = TotalN + cAllocLogFacts;
	if(m_NumLogFacts <= 2)
		{
		m_pLogFact[0] = log(1.0);
		m_pLogFact[1] = log(1.0);
		m_NumLogFacts = 2;
		}

	pLogFact = &m_pLogFact[m_NumLogFacts - 1];
	for(Idx = m_NumLogFacts; Idx < m_AllocLogFacts; Idx++,pLogFact++)
		pLogFact[1] = log((double)Idx) + *pLogFact;
	m_NumLogFacts = m_AllocLogFacts;
	}

// now to calc the one tailed P value
Num = (m_pLogFact[R1T]  + m_pLogFact[R2T] + m_pLogFact[C1T] + m_pLogFact[C2T]) - m_pLogFact[TotalN]; 
PValue = 0.0;
while(R1C1 >= 0)
	{
	PValue += exp(Num - (m_pLogFact[R1C2] + m_pLogFact[R1C1] + m_pLogFact[R2C1] + m_pLogFact[R2C2]));
	R1C1--;
	R1C2++;
	R2C1++;
	R2C2--;
	}

// with large counts then log factorials may have inherient errors so clamp returned P values to be between 0.0 and 1.0
if(PValue < 0.0)
	PValue = 0.0;
else
	if(PValue > 1.0)
		PValue = 1.0;
return(PValue);
}

#ifdef USETHISCODE
// Wald–Wolfowitz runs test
// Given a sequence containing 0's or 1's determines the number of runs
// and returns the probability of these being random
int
NumbRuns(int SeqLen,	// length of sequence to process for runs
		UINT8 *pSeq)	// sequence containing one or more runs of 0 and 1, any other values will be skipped
{
CStats Stats;
int SeqIdx;
int NumRuns;
int NumN0s;
int NumN1s;

int CurRunLen;
int MaxRunLen;

// skip any values > 1
while(SeqLen && *pSeq & 0xfe)
	{
	pSeq += 1;
	SeqLen-=1;
	}
if(!SeqLen)				
	return(-1.0);		// error

NumRuns = 1;			// initial value starts a run			 
CurRunLen = 1;			 
MaxRunLen = 1;			
NumN0s = 0;
NumN1s = 0;

for(SeqIdx = 1; SeqIdx < SeqLen; SeqIdx++, pSeq++)
	{
	if(*pSeq & 0xfe)
		continue;

	if(*pSeq == 0)
		NumN0s += 1;
	else
		NumN1s += 1;

	if(*pSeq != pSeq[1]) // terminate current and starting new run?
		{
		CurRunLen += 1;
		if(CurRunLen > MaxRunLen)
			MaxRunLen = CurRunLen;
		NumRuns++;
		CurRunLen = 1;
		continue;
		}
    CurRunLen += 1;
	}

// MaxRunLen, NumRuns, NumN0s and NumN1s are now known
int TotN = NumN0s + NumN1s;
double TWON1N2 = 2.0*NumN0s*NumN1s;
double ERuns = 1 +((TWON1N2)/TotN);	// given the number of observed N1s and N2s then ERuns is the expected number of runs
double varianceNumerator = TWON1N2 * (TWON1N2 - TotN);
double varianceDenominator = (double)((TotN*TotN)*(TotN-1));
double variance = varianceNumerator / varianceDenominator;
double Z = ((double)NumRuns - ERuns) / sqrt(variance);

// assume an exact P(r less than or equal to R)
double PrLessEqlR;
int z;
double sum = 0.0;

if(!(R&0x01))	// is even? then
	{		
	for(z = 2; z <= R; z++)
		sum += Stats.Calc_nCk(NumN0s-1,z/2-1)* Stats.Calc_nCk(NumN1s-1,z/2-1);
	PrLessEqlR = (2/Stats.Calc_nCk(NumN0s+NumN1s,NumN0s)) * sum;
	}
else		// else R is odd
	{
	for(z = 2; z <= R; z++)
		sum += (Stats.Calc_nCk(NumN0s-1,z/2-1)* Stats.Calc_nCk(NumN1s-1,z/2-1))  + Stats.Calc_nCk(NumN0s-1,k-2)* Stats.Calc_nCk(NumN1s-1,k-2));
	PrLessEqlR = (1/Stats.Calc_nCk(NumN0s+NumN1s,NumN0s)) * sum;
	}

return(PrLessEqlR);
}
#endif

// 

/* Usage: binomial(n,k,p)
 *
 * ----------------------------------------------------------------------------------------
 *
 * For n independent trials each with probability p of success, and q = 1 - p of failure,
 * the probability of k successes is given by
 *
 *	Pr(K = k) =	nCk * p^k * q^(n-k)			(for integer k)
 *
 * 		     	     n!
 *	Where nCk =	-----------  (the number of ways to select k successes from n trials)
 *			 (n-k)!k!
 *
 * The cumulative distribution is the sum of all the probabilities K = 0 up to K = k
 *
 * 			  k
 * 			 ---
 *	Pr(K <= k) =	 \	nCj * p^j * q^(n-j)		(for integer k)
 *			 /
 *			 ---
 *			j = 0
 */


// Calculates nCk = n! / (n-k)!k!
double
CStats::Calc_nCk(UINT32 n, UINT32 k) 
{
long double accum;
UINT32 Idx;
UINT64 nFactorial;
UINT64 kFactorial;
UINT64 nkFactorial;

if (k > n)
	return(0.0);

if (k > n/2)
   k = n-k; // Take advantage of symmetry

if(n < 1000)
	{
	nFactorial = 1;			
	kFactorial = 1;
	nkFactorial = 1;
	for(Idx = 2; Idx <= n; Idx++)
		{
		nFactorial *= Idx;
		if(Idx <= k)
			kFactorial *= Idx;
		if(Idx <= (n-k))
			nkFactorial *= Idx;
		}
	accum = (long double)nFactorial/(long double)(nkFactorial * kFactorial);
	}

accum = 1;
for (Idx = 1; Idx <= k; Idx++)
  accum = accum * (n-k+Idx) / Idx;

return((double)accum); 
}

/**
 * Calculates Pr(K = k) = nCk * p^k * q^(n-k)
 */
double 
CStats::ProbKeqlk(UINT32 n, UINT32 k, double p)
{
if(p < 0 || p > 1)
	return -1;
double nCk = Calc_nCk(n, k);
double p2pow = pow(p, (INT32)k);
double q2pow = pow(1 - p, (INT32)(n - k));
return (nCk * p2pow * q2pow);
}

// The cumulative distribution is the sum of all the probabilities K = 0 up to K = k
double
CStats::Binomial(int n, int k, double p)
{
if(k > n)
	return(0.0);
if(n > 5000)	// clamp n to be at most 5000 for performance and overflow reasons
	{
	k = (int)((1000.0/n) * k);
	n = 5000;
	}

//Sum probabilities over desired range
double sum = 0;
for(int i = 0; i <= k; i++)
	{
	sum += ProbKeqlk(n, i, p);
	if(sum >= 1.0)
		break;
	}
return(min(sum,1.0));
}



/* ------------------------------------------------------------------------- 
 * This is an ANSI C library that can be used to evaluate the probability 
 * density functions (pdf's), cumulative distribution functions (cdf's), and 
 * inverse distribution functions (idf's) for a variety of discrete and 
 * continuous random variables.
 *
 * The following notational conventions are used
 *                 x : possible value of the random variable
 *                 u : real variable (probability) between 0.0 and 1.0 
 *  a, b, n, p, m, s : distribution-specific parameters
 *
 * There are pdf's, cdf's and idf's for 6 discrete random variables
 *
 *      Random Variable    Range (x)  Mean         Variance
 *
 *      Bernoulli(p)       0..1       p            p*(1-p)
 *      Binomial(n, p)     0..n       n*p          n*p*(1-p)
 *      Equilikely(a, b)   a..b       (a+b)/2      ((b-a+1)*(b-a+1)-1)/12 
 *      Geometric(p)       0...       p/(1-p)      p/((1-p)*(1-p))
 *      Pascal(n, p)       0...       n*p/(1-p)    n*p/((1-p)*(1-p))
 *      Poisson(m)         0...       m            m
 *
 * and for 7 continuous random variables
 *
 *      Uniform(a, b)      a < x < b  (a+b)/2      (b-a)*(b-a)/12
 *      Exponential(m)     x > 0      m            m*m
 *      Erlang(n, b)       x > 0      n*b          n*b*b
 *      Normal(m, s)       all x      m            s*s
 *      Lognormal(a, b)    x > 0         see below
 *      Chisquare(n)       x > 0      n            2*n
 *      Student(n)         all x      0  (n > 1)   n/(n-2)   (n > 2)
 *
 * For the Lognormal(a, b), the mean and variance are
 *
 *                        mean = Exp(a + 0.5*b*b)
 *                    variance = (Exp(b*b) - 1)*Exp(2*a + b*b)
 *
 * Name            : rvms.c (Random Variable ModelS)
 * Author          : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 11-22-97
 * ------------------------------------------------------------------------- 
 */

#include <math.h>
//#include "rvms.h"

#define TINY    1.0e-10
#define SQRT2PI 2.506628274631               /* sqrt(2 * pi) */

static double pdfStandard(double x);
static double cdfStandard(double x);
static double idfStandard(double u);
static double LogGamma(double a);
static double LogBeta(double a, double b);
static double InGamma(double a, double b);
static double InBeta(double a, double b, double x);


   double pdfBernoulli(double p, long x)
/* =======================================
 * NOTE: use 0.0 < p < 1.0 and 0 <= x <= 1
 * =======================================
 */
{
   return ((x == 0) ? 1.0 - p : p);
}

   double cdfBernoulli(double p, long x)
/* =======================================
 * NOTE: use 0.0 < p < 1.0 and 0 <= x <= 1 
 * =======================================
 */
{
   return ((x == 0) ? 1.0 - p : 1.0);
}

   long idfBernoulli(double p, double u)
/* =========================================
 * NOTE: use 0.0 < p < 1.0 and 0.0 < u < 1.0 
 * =========================================
 */
{
   return ((u < 1.0 - p) ? 0 : 1);
}

   double pdfEquilikely(long a, long b, long x)
/* ============================================ 
 * NOTE: use a <= x <= b 
 * ============================================
 */
{
   return (1.0 / (b - a + 1.0));
}

   double cdfEquilikely(long a, long b, long x)
/* ============================================
 * NOTE: use a <= x <= b 
 * ============================================
 */
{
   return ((x - a + 1.0) / (b - a + 1.0));
}

   long idfEquilikely(long a, long b, double u)
/* ============================================ 
 * NOTE: use a <= b and 0.0 < u < 1.0 
 * ============================================
 */
{
   return (a + (long) (u * (b - a + 1)));
}

   double pdfBinomial(long n, double p, long x)
/* ============================================ 
 * NOTE: use 0 <= x <= n and 0.0 < p < 1.0 
 * ============================================
 */
{
   double s, t;

   s = LogChoose(n, x);
   t = x * log(p) + (n - x) * log(1.0 - p);
   return (exp(s + t));
}

   double cdfBinomial(long n, double p, long x)
/* ============================================ 
 * NOTE: use 0 <= x <= n and 0.0 < p < 1.0 
 * ============================================
 */
{
   if (x < n)
     return (1.0 - InBeta(x + 1, n - x, p));
   else
     return (1.0);
}

   long idfBinomial(long n, double p, double u)
/* ================================================= 
 * NOTE: use 0 <= n, 0.0 < p < 1.0 and 0.0 < u < 1.0 
 * =================================================
 */
{
   long x = (long) (n * p);             /* start searching at the mean */

   if (cdfBinomial(n, p, x) <= u)
     while (cdfBinomial(n, p, x) <= u)
       x++;
   else if (cdfBinomial(n, p, 0) <= u)
     while (cdfBinomial(n, p, x - 1) > u)
       x--;
   else
     x = 0;
   return (x);
}

   double pdfGeometric(double p, long x)
/* ===================================== 
 * NOTE: use 0.0 < p < 1.0 and x >= 0 
 * =====================================
 */
{
   return ((1.0 - p) * exp(x * log(p)));
}

   double cdfGeometric(double p, long x)
/* ===================================== 
 * NOTE: use 0.0 < p < 1.0 and x >= 0 
 * =====================================
 */
{
   return (1.0 - exp((x + 1) * log(p)));
}

   long idfGeometric(double p, double u)
/* ========================================= 
 * NOTE: use 0.0 < p < 1.0 and 0.0 < u < 1.0 
 * =========================================
 */
{
   return ((long) (log(1.0 - u) / log(p)));
}

   double pdfPascal(long n, double p, long x)
/* =========================================== 
 * NOTE: use n >= 1, 0.0 < p < 1.0, and x >= 0 
 * ===========================================
 */
{
   double  s, t;

   s = LogChoose(n + x - 1, x);
   t = x * log(p) + n * log(1.0 - p);
   return (exp(s + t));
}

   double cdfPascal(long n, double p, long x)
/* =========================================== 
 * NOTE: use n >= 1, 0.0 < p < 1.0, and x >= 0 
 * ===========================================
 */
{
   return (1.0 - InBeta(x + 1, n, p));
}

   long idfPascal(long n, double p, double u)
/* ================================================== 
 * NOTE: use n >= 1, 0.0 < p < 1.0, and 0.0 < u < 1.0 
 * ==================================================
 */
{
   long x = (long) (n * p / (1.0 - p));    /* start searching at the mean */

   if (cdfPascal(n, p, x) <= u)
     while (cdfPascal(n, p, x) <= u)
       x++;
   else if (cdfPascal(n, p, 0) <= u)
     while (cdfPascal(n, p, x - 1) > u)
       x--;
   else
     x = 0;
   return (x);
}

   double pdfPoisson(double m, long x)
/* ===================================
 * NOTE: use m > 0 and x >= 0 
 * ===================================
 */
{
   double t;

   t = - m + x * log(m) - LogFactorial(x);
   return (exp(t));
}

   double cdfPoisson(double m, long x)
/* =================================== 
 * NOTE: use m > 0 and x >= 0 
 * ===================================
 */
{
   return (1.0 - InGamma(x + 1, m));
}

   long idfPoisson(double m, double u)
/* =================================== 
 * NOTE: use m > 0 and 0.0 < u < 1.0 
 * ===================================
 */
{
   long x = (long) m;                    /* start searching at the mean */

   if (cdfPoisson(m, x) <= u)
     while (cdfPoisson(m, x) <= u)
       x++;
   else if (cdfPoisson(m, 0) <= u)
     while (cdfPoisson(m, x - 1) > u)
       x--;
   else
     x = 0;
   return (x);
}

   double pdfUniform(double a, double b, double x)
/* =============================================== 
 * NOTE: use a < x < b 
 * ===============================================
 */
{
   return (1.0 / (b - a));
}

   double cdfUniform(double a, double b, double x)
/* =============================================== 
 * NOTE: use a < x < b 
 * ===============================================
 */
{
   return ((x - a) / (b - a));
}

   double idfUniform(double a, double b, double u)
/* =============================================== 
 * NOTE: use a < b and 0.0 < u < 1.0 
 * ===============================================
 */
{
   return (a + (b - a) * u);
}

   double pdfExponential(double m, double x)
/* ========================================= 
 * NOTE: use m > 0 and x > 0 
 * =========================================
 */
{
   return ((1.0 / m) * exp(- x / m));
}

   double cdfExponential(double m, double x)
/* ========================================= 
 * NOTE: use m > 0 and x > 0 
 * =========================================
 */
{
   return (1.0 - exp(- x / m));
}

   double idfExponential(double m, double u)
/* ========================================= 
 * NOTE: use m > 0 and 0.0 < u < 1.0 
 * =========================================
 */
{
   return (- m * log(1.0 - u));
}

   double pdfErlang(long n, double b, double x)
/* ============================================ 
 * NOTE: use n >= 1, b > 0, and x > 0 
 * ============================================
 */
{
   double t;

   t = (n - 1) * log(x / b) - (x / b) - log(b) - LogGamma(n);
   return (exp(t));
}

   double cdfErlang(long n, double b, double x)
/* ============================================ 
 * NOTE: use n >= 1, b > 0, and x > 0 
 * ============================================
 */
{
   return (InGamma(n, x / b));
}

   double idfErlang(long n, double b, double u)
/* ============================================ 
 * NOTE: use n >= 1, b > 0 and 0.0 < u < 1.0 
 * ============================================
 */
{
   double t, x = n * b;                   /* initialize to the mean, then */

   do {                                   /* use Newton-Raphson iteration */
     t = x;
     x = t + (u - cdfErlang(n, b, t)) / pdfErlang(n, b, t);
     if (x <= 0.0)
       x = 0.5 * t;
   } while (fabs(x - t) >= TINY);
   return (x);
}

   static double pdfStandard(double x)
/* =================================== 
 * NOTE: x can be any value 
 * ===================================
 */
{
   return (exp(- 0.5 * x * x) / SQRT2PI);
}

   static double cdfStandard(double x)
/* =================================== 
 * NOTE: x can be any value 
 * ===================================
 */
{ 
   double t;

   t = InGamma(0.5, 0.5 * x * x);
   if (x < 0.0)
     return (0.5 * (1.0 - t));
   else
     return (0.5 * (1.0 + t));
}

   static double idfStandard(double u)
/* =================================== 
 * NOTE: 0.0 < u < 1.0 
 * ===================================
 */
{ 
   double t, x = 0.0;                    /* initialize to the mean, then  */

   do {                                  /* use Newton-Raphson iteration  */
     t = x;
     x = t + (u - cdfStandard(t)) / pdfStandard(t);
   } while (fabs(x - t) >= TINY);
   return (x);
}

   double pdfNormal(double m, double s, double x)
/* ============================================== 
 * NOTE: x and m can be any value, but s > 0.0 
 * ==============================================
 */
{ 
   double t = (x - m) / s;

   return (pdfStandard(t) / s);
}

   double cdfNormal(double m, double s, double x)
/* ============================================== 
 * NOTE: x and m can be any value, but s > 0.0 
 * ==============================================
 */
{ 
   double t = (x - m) / s;

   return (cdfStandard(t));
}

   double idfNormal(double m, double s, double u)
/* ======================================================= 
 * NOTE: m can be any value, but s > 0.0 and 0.0 < u < 1.0 
 * =======================================================
 */
{
   return (m + s * idfStandard(u));
}

   double pdfLognormal(double a, double b, double x)
/* =================================================== 
 * NOTE: a can have any value, but b > 0.0 and x > 0.0 
 * ===================================================
 */
{ 
   double t = (log(x) - a) / b;

   return (pdfStandard(t) / (b * x));
}

   double cdfLognormal(double a, double b, double x)
/* =================================================== 
 * NOTE: a can have any value, but b > 0.0 and x > 0.0 
 * ===================================================
 */
{ 
   double t = (log(x) - a) / b;

   return (cdfStandard(t));
}

   double idfLognormal(double a, double b, double u)
/* ========================================================= 
 * NOTE: a can have any value, but b > 0.0 and 0.0 < u < 1.0 
 * =========================================================
 */
{ 
   double t;

   t = a + b * idfStandard(u);
   return (exp(t));
}

   double pdfChisquare(long n, double x)
/* ===================================== 
 * NOTE: use n >= 1 and x > 0.0 
 * =====================================
 */
{ 
   double t, s = n / 2.0;

   t = (s - 1.0) * log(x / 2.0) - (x / 2.0) - log(2.0) - LogGamma(s);
   return (exp(t));
}

   double cdfChisquare(long n, double x)
/* ===================================== 
 * NOTE: use n >= 1 and x > 0.0 
 * =====================================
 */
{
   return (InGamma(n / 2.0, x / 2));
}

   double idfChisquare(long n, double u)
/* ===================================== 
 * NOTE: use n >= 1 and 0.0 < u < 1.0 
 * =====================================
 */
{ 
   double t, x = n;                         /* initialize to the mean, then */

   do {                                     /* use Newton-Raphson iteration */
     t = x;
     x = t + (u - cdfChisquare(n, t)) / pdfChisquare(n, t);
     if (x <= 0.0)
       x = 0.5 * t;
   } while (fabs(x - t) >= TINY);
   return (x);
}

   double pdfStudent(long n, double x)
/* =================================== 
 * NOTE: use n >= 1 and x > 0.0 
 * ===================================
 */
{ 
   double s, t;

   s = -0.5 * (n + 1) * log(1.0 + ((x * x) / (double) n));
   t = -LogBeta(0.5, n / 2.0);
   return (exp(s + t) / sqrt((double) n));
}

   double cdfStudent(long n, double x)
/* =================================== 
 * NOTE: use n >= 1 and x > 0.0 
 * ===================================
 */
{ 
   double s, t;

   t = (x * x) / (n + x * x);
   s = InBeta(0.5, n / 2.0, t);
   if (x >= 0.0)
     return (0.5 * (1.0 + s));
   else
     return (0.5 * (1.0 - s));
}

   double idfStudent(long n, double u)
/* =================================== 
 * NOTE: use n >= 1 and 0.0 < u < 1.0 
 * ===================================
 */
{ 
   double t, x = 0.0;                       /* initialize to the mean, then */

   do {                                     /* use Newton-Raphson iteration */
     t = x;
     x = t + (u - cdfStudent(n, t)) / pdfStudent(n, t);
   } while (fabs(x - t) >= TINY);
   return (x);
}

/* ===================================================================
 * The six functions that follow are a 'special function' mini-library
 * used to support the evaluation of pdf, cdf and idf functions.
 * ===================================================================
 */

   static double LogGamma(double a)
/* ======================================================================== 
 * LogGamma returns the natural log of the gamma function.
 * NOTE: use a > 0.0 
 *
 * The algorithm used to evaluate the natural log of the gamma function is
 * based on an approximation by C. Lanczos, SIAM J. Numerical Analysis, B,
 * vol 1, 1964.  The constants have been selected to yield a relative error
 * which is less than 2.0e-10 for all positive values of the parameter a.    
 * ======================================================================== 
 */
{ 
   double s[6], sum, temp;
   int    i;

   s[0] =  76.180091729406 / a;
   s[1] = -86.505320327112 / (a + 1.0);
   s[2] =  24.014098222230 / (a + 2.0);
   s[3] =  -1.231739516140 / (a + 3.0);
   s[4] =   0.001208580030 / (a + 4.0);
   s[5] =  -0.000005363820 / (a + 5.0);
   sum  =   1.000000000178;
   for (i = 0; i < 6; i++) 
     sum += s[i];
   temp = (a - 0.5) * log(a + 4.5) - (a + 4.5) + log(SQRT2PI * sum);
   return (temp);
}

   double LogFactorial(long n)
/* ==================================================================
 * LogFactorial(n) returns the natural log of n!
 * NOTE: use n >= 0
 *
 * The algorithm used to evaluate the natural log of n! is based on a
 * simple equation which relates the gamma and factorial functions.
 * ==================================================================
 */
{
   return (LogGamma(n + 1));
}

   static double LogBeta(double a, double b)
/* ======================================================================
 * LogBeta returns the natural log of the beta function.
 * NOTE: use a > 0.0 and b > 0.0
 *
 * The algorithm used to evaluate the natural log of the beta function is 
 * based on a simple equation which relates the gamma and beta functions.
 *
 */
{ 
   return (LogGamma(a) + LogGamma(b) - LogGamma(a + b));
}

   double LogChoose(long n, long m)
/* ========================================================================
 * LogChoose returns the natural log of the binomial coefficient C(n,m).
 * NOTE: use 0 <= m <= n
 *
 * The algorithm used to evaluate the natural log of a binomial coefficient
 * is based on a simple equation which relates the beta function to a
 * binomial coefficient.
 * ========================================================================
 */
{
   if (m > 0)
     return (-LogBeta(m, n - m + 1) - log((double)m));
   else
     return (0.0);
}

   static double InGamma(double a, double x)
/* ========================================================================
 * Evaluates the incomplete gamma function.
 * NOTE: use a > 0.0 and x >= 0.0
 *
 * The algorithm used to evaluate the incomplete gamma function is based on
 * Algorithm AS 32, J. Applied Statistics, 1970, by G. P. Bhattacharjee.
 * See also equations 6.5.29 and 6.5.31 in the Handbook of Mathematical
 * Functions, Abramowitz and Stegum (editors).  The absolute error is less 
 * than 1e-10 for all non-negative values of x.
 * ========================================================================
 */
{ 
   double t, sum, term, factor, f, g, c[2], p[3], q[3];
   long   n;

   if (x > 0.0)
     factor = exp(-x + a * log(x) - LogGamma(a));
   else
     factor = 0.0;
   if (x < a + 1.0) {                 /* evaluate as an infinite series - */
     t    = a;                        /* A & S equation 6.5.29            */
     term = 1.0 / a;
     sum  = term;
     while (term >= TINY * sum) {     /* sum until 'term' is small */
       t++;
       term *= x / t;
       sum  += term;
     } 
     return (factor * sum);
   }
   else {                             /* evaluate as a continued fraction - */
     p[0]  = 0.0;                     /* A & S eqn 6.5.31 with the extended */
     q[0]  = 1.0;                     /* pattern 2-a, 2, 3-a, 3, 4-a, 4,... */
     p[1]  = 1.0;                     /* - see also A & S sec 3.10, eqn (3) */
     q[1]  = x;
     f     = p[1] / q[1];
     n     = 0;
     do {                             /* recursively generate the continued */
       g  = f;                        /* fraction 'f' until two consecutive */
       n++;                           /* values are small                   */
       if ((n % 2) > 0) {
         c[0] = ((double) (n + 1) / 2) - a;
         c[1] = 1.0;
       }
       else {
         c[0] = (double) n / 2;
         c[1] = x;
       }
       p[2] = c[1] * p[1] + c[0] * p[0];
       q[2] = c[1] * q[1] + c[0] * q[0];
       if (q[2] != 0.0) {             /* rescale to avoid overflow */
         p[0] = p[1] / q[2];
         q[0] = q[1] / q[2];
         p[1] = p[2] / q[2];
         q[1] = 1.0;
         f    = p[1];
       }
     } while ((fabs(f - g) >= TINY) || (q[1] != 1.0));
     return (1.0 - factor * f);
   }
}

   static double InBeta(double a, double b, double x)
/* ======================================================================= 
 * Evaluates the incomplete beta function.
 * NOTE: use a > 0.0, b > 0.0 and 0.0 <= x <= 1.0
 *
 * The algorithm used to evaluate the incomplete beta function is based on
 * equation 26.5.8 in the Handbook of Mathematical Functions, Abramowitz
 * and Stegum (editors).  The absolute error is less than 1e-10 for all x
 * between 0 and 1.
 * =======================================================================
 */
{ 
   double t, factor, f, g, c, p[3], q[3];
   int    swap;
   long   n;

   if (x > (a + 1.0) / (a + b + 1.0)) { /* to accelerate convergence   */
     swap = 1;                          /* complement x and swap a & b */
     x    = 1.0 - x;
     t    = a;
     a    = b;
     b    = t;
   }
   else                                 /* do nothing */
     swap = 0;
   if (x > 0)
     factor = exp(a * log(x) + b * log(1.0 - x) - LogBeta(a,b)) / a;
   else
     factor = 0.0;
   p[0] = 0.0;
   q[0] = 1.0;
   p[1] = 1.0;
   q[1] = 1.0;
   f    = p[1] / q[1];
   n    = 0;
   do {                               /* recursively generate the continued */
     g = f;                           /* fraction 'f' until two consecutive */
     n++;                             /* values are small                   */
     if ((n % 2) > 0) {
       t = (double) (n - 1) / 2;
       c = -(a + t) * (a + b + t) * x / ((a + n - 1.0) * (a + n));
     }
     else {
       t = (double) n / 2;
       c = t * (b - t) * x / ((a + n - 1.0) * (a + n));
     }
     p[2] = p[1] + c * p[0];
     q[2] = q[1] + c * q[0];
     if (q[2] != 0.0) {                 /* rescale to avoid overflow */
       p[0] = p[1] / q[2];
       q[0] = q[1] / q[2];
       p[1] = p[2] / q[2];
       q[1] = 1.0;
       f    = p[1];
     }
   } while ((fabs(f - g) >= TINY) || (q[1] != 1.0));
   if (swap) 
     return (1.0 - factor * f);
   else
     return (factor * f);
}

