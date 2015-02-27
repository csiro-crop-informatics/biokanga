#pragma once

const int cBSFRdsVersion = 5;			// latest file header version which can be handled
const int cBSFRdsVersionBack= 3;		// backward compatible to this version

const int cMaxWorkerThreads = 128;		// limiting max number of threads to this many
const int cMaxSubseqsPerBlock = 0x07fff;// max number of subsequences for processing per thread as a block

const int cDfltAllowedSubs = 2;			// default number of aligner induced substitutions
const int cMaxAllowedSubs = 30;			// allow at most this many aligner induced substitutions

const int cDfltMaxMatches = 5000;		// default : allow at most this many subsequence matches before filtering out the subsequence
const int cMaxMaxMatches = 1000000;		// allow at most this many subsequence matches before filtering out the subsequence

const int cMinCoreLen = 5;				// absolute minimum core length supported
const int cMaxCoreLen = 100;			// maximum core length supported
const double cDfltzygosityThreshold = 0.25f; // default Zygosity threshold
const int cDfltSensCoreIters  = 5000;	// default sensitivity core explore depth 
const int cMoreSensCoreIters  = 25000;	// more sensitivity core explore depth
const int cUltraSensCoreIters = 50000;	// ultra sensitivity core explore depth
const int cMinSensCoreIters   = 1000;	// min sensitivity core explore depth

const int cDfltMaxNs = 1;				// default max number of indeterminate 'N's in the read or aligned to subsequence
const int cMaxNs = 5;					// allow user to specify at most this number of 'N's in the read or aligned to subsequence

const int cMinSubseqLen = 20;			// allow for subsequence sample lengths >= cMinSubseqLen
const int cDlftSubseqLen = 25;			// default subsequence sample length == cDlftSubseqLen
const int cMaxSubseqLen = 1000;			// allow for subsequence sample length <= cMaxSubseqLen			

const int cAllocLineBuffSize = 8196000; // when writing to results file then can buffer upto this many chars

