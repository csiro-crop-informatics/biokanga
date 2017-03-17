Version 4.3.x
===========================
Released March 1st 2017

BioKanga is an integrated toolkit of high performance bioinformatics subprocesses targeting the challenges of next generation sequencing analytics. Kanga is an acronym standing for 'K-mer Adaptive Next Generation Aligner'.

Why YAL (Yet Another Aligner)
-----------------------------
The BioKanga alignment subprocess is a highly efficient short-read aligner which incorporates an empirically derived understanding of sequence uniqueness within a target genome to enable robust alignment of next generation sequencer short read datasets in either colorspace (ABI SOLiD) or basespace (Illumina).
Compared with other widely used aligners, BioKanga provides substantial gains in both the proportion and quality of aligned sequence reads at competitive or increased computational efficiency. Unlike most other aligners, BioKanga utilises Hamming distances between putative alignments to the targeted genome assembly for any given read as the discrimative acceptance criteria rather than relying on sequencer generated quality scores.

Another primary differentiator for BioKanga is that this toolkit can process billions of reads against targeted genomes containing 100 million contigs and totalling up to 100Gbp of sequence.

Availability
------------
Binary pre-built releases are available for either Linux or Windows x64 hosting environments, and source code (C++) is also available for those requiring a local build.

Toolset Components
------------------
The BioKanga toolset contains a number of subprocesses, each of which is targeting a specific bioinformatics analytics task. Primary subprocesses provide functionality for:
 - Generate simulated NGS datasets
 - Quality check the raw NGS reads to identify potential processing issues
 - Filter NGS reads for sequencer errors and/or exact duplicates
 - de Novo assemble filtered reads into contigs
 - Scaffold de Novo assembled contigs
 - Blitz local alignments
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
 - Remap alignment loci
 - Locate and report regions of interest
 - Generate marker sequences from SNP loci
 - Generate SQLite Marker Database from SNP markers
 - Generate SQLite SNP Database from aligner identified SNPs
 - Generate SQLite DE Database from RNA-seq DE
 - Generate SQLite Blat alignment PSL database

Installation
------------
To install, simply copy the BioKanga process binary image (biokanga) into a directory which is on your executable path. There are a number of user specified parameters specific to each subprocess of the BioKanga toolset. Most of these are optional, and optional parameters would generally only be required to be specified by the user when targeting some specific alignment issue. In general, the only mandatory parameters are those specifying input datasets and where to write output result sets.

Developers
----------
BioKanga is actively being developed and enhanced by Dr Stuart Stephen, with contributions from other group members, in Dr Alex Whan's Bioinformatics Team at the CSIRO, Canberra, Australia.

Please report issues and comments to: alex.whan@csiro.au


























