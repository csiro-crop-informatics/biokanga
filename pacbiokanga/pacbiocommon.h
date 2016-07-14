#pragma once
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


const int cMinSeedCoreLen = 12;							// user can specify seed cores down to this minimum length
const int cDfltSeedCoreLen = 14;						// default seed cores of this length
const int cMaxSeedCoreLen = 50;							// user can specify seed cores up to to this maximum length
const int cDfltScaffSeedCoreLen = 35;					// default seed cores of this length when generating overlap detail for scaffolding
const int cDfltConsolidateSeedCoreLen = 25;             // default seed cores of this length when consolidating transcripts into a reference transcript

const int cDfltDeltaCoreOfs = 2;						// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps
const int cDfltConsolidateDeltaCoreOfs = 5;				// offset by this many bp the core windows of coreSeqLen along the probe sequence when checking for overlaps if consolidating transcripts
const int cMaxDeltaCoreOfs = 25;						// allowing at most this delta core offset

const int cDfltMaxSeedCoreDepth = 15000;                // only extend a seed core if there are no more than this number of matching cores in all targeted sequences
const int cDfltMaxConsolidateSeedCoreDepth = 75000;     // only extend a seed core if there are no more than this number of matching cores in all targeted sequences if consolidating transcripts

const int cDfltMaxAcceptHitsPerSeedCore = 10000;		// only accept up to this many extended cores at any probe offset 
const int cMinNumSeedCores = 1;							// user can specify requiring at least this many seed cores between overlapping scaffold sequences
const int cDfltNumSeedCores = 10;						// default is to require at least this many seed cores between overlapping scaffold sequences
const int cDfltConsolidateNumSeedCores = 15;			// default is to require at least this many seed cores between overlapping sequences when consolidating transcripts into a reference transcript
const int cDfltScaffSeedCores = 30;				        // default is to require at least this many seed cores between overlapping scaffold sequences if generating overlap details, sequences expected to be relatively error free 
const int cMaxNumSeedCores = 50;						// user can specify requiring up to many seed cores between overlapping scaffold sequences
const int cAnchorLen = 6;								// require 5' and 3' end anchors of at least this length for overlap merging
const int cQualCoreKMerLen = 3;							// using tri-mers when checking for core downstream shared kmers between probe and target. Note currently a max of 4 would be supported as any more would violate cQualCoreDelta constraint 
const int cQualCoreDelta = (cQualCoreKMerLen * 2) + 1;	 // looking for matching trimers starting within +/- 7bp of the probe trimer. NOTE must be <= cMinSeedCoreLen
const int cQualCoreThres = 25;							 // require at least this many kmers to be shared between probe and target before accepting core
const int cQualCoreHomopolymer = 90;                     // if any core contains more than this percentage of the same base then treat as being a near homopolymer core (likely a PacBio insert) and slough

const int cMinPBSeqLen = 500;						 // allowing for minimum PacBio sequences to be specified down to this length (could be targeting RNA transcriptome)
const int cDfltMinPBSeqLen = 10000;					 // default is to allow for error correction of PacBio sequences down to this minimum length
const int cDfltMinConsolidatePBSeqLen = 750;		 // default is to allow for PacBio sequences down to this minimum length if consolidating
const int cDfltMaxPBSeqLen = 35000;					 // default is to allow for PacBio sequences up to to this maximum length
const int cDfltMinErrCorrectLen = 5000;              // default is for this minimum length error corrected PacBio sequences
const int cMaxMinPBSeqLen = 100000;					 // allowing for minimum PacBio sequences to be specified up to this length 
const int cMaxMaxPBSeqLen = 500000;					 // allowing for maximum PacBio sequences to be specified up to this length 

const int cDfltSWMatchScore = 3;						// default SW match score for pacbio alignments
const int cDfltConsolidateSWMatchScore = 1;				// default SW match score for pacbio alignments if consolidating

const int cDfltSWMismatchPenalty = -7;					// default SW mismatch penalty for pacbio alignments
const int cDfltConsolidateSWMismatchPenalty = -10;		// default SW mismatch penalty for pacbio alignments if consolidating

const int cDfltSWGapOpenPenalty = -4;					// default SW gap open penalty for pacbio alignments
const int cDfltConsolidateSWGapOpenPenalty = -8;		// default SW gap open penalty for pacbio alignments if consolidating

const int cDfltSWGapExtnPenalty = -1;					// default SW gap extension penalty for pacbio alignments
const int cDfltConsolidateSWGapExtnPenalty = -3;		// default SW gap extension penalty for pacbio alignments if consolidating

const int cDfltSWProgExtnLen = 2;						// default SW gap extension penalties apply with gaps of at least this size
const int cDfltConsolidateSWProgExtnLen = 2;			// default SW gap extension penalties apply with gaps of at least this size if consolidating

const int cMaxAllowedSWScore = 50;						// allow SW scores or penalties to be specified up to this max

const int cMinSWPeakScore = 50;							// SW alignments must have peak scores of at least this
const int cMinSWAlignLen = 50;							// and alignment length of at least this many bases to be further processed

const int cDfltMaxOverlapFloat = 1500;					// allow up to this much float on overlaps to account for the PacBio error profile 
const int cDfltScaffMaxOverlapFloat = 200;				// but when assembling or correcting contigs with error corrected reads then reduce
const int cDfltConsolidateMaxOverlapFloat = 100;		// allow up to this much float on overlaps to account for the PacBio error profile when consolidating 


const int cAllocdNumCoreHits = 1000000;					 // each thread preallocs for this many core hits, realloc'd as may be required
const int cAllocdQuerySeqLen = 500000;					 // each thread preallocs to hold query sequences of this length, realloc'd as may be required
const int cSummaryTargCoreHitCnts = 20000;				 // summary core hit counts on at most this many targets

const int cChkOverlapGapLen = 20;						 // if gap between probe cores with at least one match > this threshold then set core probe offset at which to check for core extension (set 0 to disable core extensions)

const int cMaxWorkerThreads = 128;							// limiting max number of threads to this many

const UINT32 cMaxValidID = 0xffffff00;      // treat any vertex, edge, subgraph or sequence identifier over this value as a processing error identifier which is to be cast as teBSFrsltCodes for actual error 
typedef UINT32 tVertID;						// vertice identifiers are 32bits with identifiers > cMaxValidID used as processing error indicator, cast to teBSFrsltCodes for actual error
typedef UINT32 tEdgeID;						// edges identifiers are 32bits with identifiers > cMaxValidID used as processing error indicator, cast to teBSFrsltCodes for actual error
typedef UINT32 tSeqID;						// sequence identifiers are 32bits with identifiers > cMaxValidID used as processing error indicator, cast to teBSFrsltCodes for actual error
typedef UINT32 tComponentID;			    // to contain disconnected subgraph identifiers with identifiers > cMaxValidID used as processing error indicator, cast to teBSFrsltCodes for actual error