Version 3.7.6
===========================
Released May 11th 2015

BioKanga is an integrated toolkit of high performance bioinformatics subprocesses
targeting the challenges of next generation sequencing analytics. Kanga is an
acronym standing for 'K-mer Adaptive Next Generation Aligner'.

The initial release (2.76.2) was a major release which integrated a number of
previously standalone processes into a single process executable and
incorporates many functional enhancements.

Release (2.95.0) added additional functionality to Biokanga:
 - new subprocess analysis of raw readsets for quality assurance and reports a
   number of characterisations including K-Mer distributions, K-mer
   concordances, duplicated read rates, and Phred quality score distributions.
 - putative contaminant ( e.g. adaptor/primer) sequences are processed for
   partial overlaps onto reads. This new functionality is incorporated into the
   quality, alignment, and filtering subprocesses.

Release 2.97.0 contains improvements primarily within the de Novo assembly
and scaffolding subprocesses. An additional subprocess allowing for simulation
of NGS readsets has been incorporated into the Biokanga toolset.

Release 2.97.1 is a bug fix for edge case when scaffolding.

Release 2.97.4 increases de Novo assembly througput by reducing the
default number of required steps over which the K-mer overlap is
progressively reduced, and increases the number of allowed sample files to 75
for control and experiment in the gendeseq subprocess.

Release 2.98.0
Limit of number of input files raised to 100 for SE and 200 (PE1 and PE2 each
100) files accepted. Rate of intermediate log progress writes reduced to 1
every 10 minutes (formerly was 1 every minute) and main thread stack size in
sub-process 'scaffold' increased to 100MB.

Release 2.98.8
Added subprocess to import Blat generated PSL alignments into SQLite, and
updated the NGSqc process to automatically generate SVG summary plots utilising
the PLPlot functionality.

Release 3.0.4
SSR detection subprocess now reports flanking sequences around identified
SSRs, and subprocess 'fasta2nxx' added to report length distributions as
N10..N90s over sequences in specified fasta files.

Release 3.0.9
Improved scaffolding and de Novo assembly processing.

Release 3.1.0
Extended the support for mate pair end orientations in de Novo assembly and
scaffolding.

Release 3.1.1
Iterative scaffolding with different mate pair libraries supported.

Release 3.3.5
Add new 'blitz' DNA seeded local alignment processing module - like
BLAT but much higher throughput.

Release 3.4.0
Integrated previously standalone processes for reporting regions of
interest which have coverage by read alignments (locateroi) and remapping of
read alignment loci (remaploci).

Release 3.4.2
Added vector contaminate processing to ngsqc and filter submodules

Release 3.4.3
Implemented auto-switch to utilising BAM CSI indexes if any targeted reference
sequence is longer than the 512Mbp limit imposed by BAM SAI indexes.

Version 3.4.5
Incremental changes - read lengths can now be up to 16Kbp, seed contigs for
assembly can be fastq (16Kbp max), and default insert size range changed to be
min 100bp and max 1000bp inclusive when aligning PE

Version 3.4.18
Incremental changes - read lengths further extended to 64Kbp, and 'biokanga
index' can generate a psueudorandom genome of up to 1Tbp and index same for
benchmarking on large memory systems

Version 3.5.1
Incremental changes - filtering further optimised for throughput performance
and assembly PE processing now treats individual ends as if SE in final
phases.

Version 3.5.9
Mainly a bug fix for the filter subprocess as since release 3.4.0 support for
Phred score filtering was inadvertently only partially supported due to lack
of user demand for the feature. This release fully restores Phred score
filtering functionality.

Version 3.6.0
Bug fix for NGSQC segfault under Ubuntu 14.04 when compiled optimised
and Biokanga now under GIT source control

Version 3.7.6
Added chimeric read alignment capability and the ability to call di/tri-SNPs
to the aligner.

Why YAL (Yet Another Aligner)
-----------------------------
The BioKanga alignment subprocess is a highly efficient short-read aligner
which incorporates an empirically derived understanding of sequence uniqueness
within a target genome to enable robust alignment of next generation sequencer
short read datasets in either colorspace (ABI SOLiD) or basespace (Illumina).
Compared with other widely used aligners, BioKanga provides substantial gains
in both the proportion and quality of aligned sequence reads at competitive or
increased computational efficiency. Unlike most other aligners, BioKanga
utilises Hamming distances between putative alignments to the targeted genome
assembly for any given read as the discrimative acceptance criteria rather
than relying on sequencer generated quality scores.

Another primary differentiator for BioKanga is that this toolkit can process
billions of reads against targeted genomes containing 100 million contigs and
totalling up to 100Gbp of sequence.

Availability
------------
Binary pre-built releases are available for either Linux or Windows x64 hosting
environments, and source code (C++) is also available for those requiring a
local build.

Toolset Components
------------------
The BioKanga toolset contains a number of subprocesses, each of which is
targeting a specific bioinformatics analytics task. Primary subprocesses
provide functionality for:
 - Generate simulated NGS datasets
 - Quality check the raw NGS reads to identify potential processing issues
 - Filter NGS reads for sequencer errors and/or exact duplicates
 - de Novo assemble filtered reads into contigs
 - Scaffold de Novo assembled contigs
 - Blitz DNA seeded local alignments - like BLAT but much higher throughput
 - Generate index over genome assembly or sequences
 - NGS reads alignment-less K-mer derived marker sequences generation
 - NGS reads alignment-less prefix K-mer derived marker sequences generation
 - Concatenate sequences to create pseudo-genome assembly
 - Align NGS reads to indexed genome assembly or sequences
 - Scaffold assembly contigs using PE read alignments
 - Identify SSRs in multifasta sequences
 - Map aligned reads loci to known features
 - RNA-seq differential expression analyser with optional Pearsons generation
 - Generate tab delimited counts file for input to DESeq or EdgeR
 - Extract fasta sequences from multifasta file
 - Merge PE short insert overlap reads
 - SNP alignment derived marker sequences identification
 - Reporting of identified regions of interest
 - Remapping of alignment loci
 - Generate marker sequences from SNP loci
 - Generate SQLite Marker Database from SNP markers
 - Generate SQLite SNP Database from aligner identified SNPs
 - Generate SQLite DE Database from RNA-seq DE
 - Generate SQLite Blat alignment PSL database

Installation
------------
To install, simply copy the BioKanga process binary image (biokanga) into a
directory which is on your executable path. There are a number of user
specified parameters specific to each subprocess of the BioKanga toolset. Most
of these are optional, and optional parameters would generally only be required
to be specified by the user when targeting some specific alignment issue.
In general, the only mandatory parameters are those specifying input datasets
and where to write output result sets.

Developers
----------
BioKanga is actively being developed and enhanced by Dr Stuart Stephen, with
contributions from other group members, in Dr Jen Taylor's Bioinformatics
Group at the CSIRO, Canberra, Australia.

Please report issues and comments to: stuart.stephen@csiro.au
