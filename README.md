[![Travis Build Status](https://travis-ci.org/csiro-crop-informatics/biokanga.svg?branch=master)](https://travis-ci.org/csiro-crop-informatics/biokanga)
[![Build status](https://ci.appveyor.com/api/projects/status/7087y5pwrb8va0uv/branch/master?svg=true)](https://ci.appveyor.com/project/alexwhan/biokanga/branch/master)
[![Docker pulls](https://img.shields.io/docker/build/csirocropinformatics/biokanga.svg?logo=docker)](https://hub.docker.com/r/csirocropinformatics/biokanga)
[![Docker pulls](https://img.shields.io/docker/automated/csirocropinformatics/biokanga.svg?logo=docker)](https://hub.docker.com/r/csirocropinformatics/biokanga)
[![Docker pulls](https://img.shields.io/docker/pulls/csirocropinformatics/biokanga.svg?logo=docker)](https://hub.docker.com/r/csirocropinformatics/biokanga)



# BioKanga 
BioKanga is an integrated toolkit of high performance bioinformatics subprocesses targeting the challenges of next generation sequencing analytics. Kanga is an acronym standing for 'K-mer Adaptive Next Generation Aligner'.

## Why YAL (Yet Another Aligner)
Compared with other widely used aligners, BioKanga provides substantial gains in both the proportion and quality of aligned sequence reads at competitive or increased computational efficiency. Unlike most other aligners, BioKanga utilises Hamming distances between putative alignments to the targeted genome assembly for any given read as the discrimative acceptance criteria rather than relying on sequencer generated quality scores.

Another primary differentiator for BioKanga is that this toolkit can process billions of reads against targeted genomes containing 100 million contigs and totalling up to 100Gbp of sequence.

## Toolset Components
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


## Build and installation
### Linux
To build on linux, clone this repository, run `autoreconf`, `configure` and `make`. The following example will install the biokanga toolkit to a `bin` directory underneath the user's home directory.

```
git clone https://github.com/csiro-crop-informatics/biokanga.git
cd biokanga
autoreconf -f -i
./configure --prefix=$HOME
make install
```

Alternatively, the binary built for the appropriate platform can be used directly.

### Windows
To build on Windows, the current version requires Visual Studio 2015 or 2017 **with build tools v140**. 
1. Open the `biokanga.sln` file in Visual Studio. 
2. Under the Build menu, select Configuration Manager. 
3. For Active solution platform, select x64. 
4. The project can then be built. By default, executables will be copied into the `Win64` directory.

Alternatively, the windows binaries can be used directly.

## Documentation
Documentation for the core functionality of biokanga and pacbiokanga is available under the `Docs` directory.

## Contributing
BioKanga is maintained by the Crop Bioinformatics and Data Science team at CSIRO in Canberra, Australia. 

Contributions are most welcome. To contribute, follow these steps.

1. Fork biokanga into your own repository ([more information](https://help.github.com/articles/about-forks/))
2. Clone and enter the repository to your development machine
3. Checkout the `dev` branch
4. Make and checkout a new branch for your work (`git checkout -b great-new-feature`)
5. Make regular commits on your new branch
6. Push your branch back to your github repository (`git push origin great-new-feature`)
7. Create a pull request to the `dev` branch of the csiro-crop-informatics/biokanga repository ([more information](https://help.github.com/articles/creating-a-pull-request/))
8. If you're work is related to an existing issue, refer to the issue in the pull request comment


## Issues
Please report issues on the [github project](https://github.com/csiro-crop-informatics/biokanga/issues).

## Authors
BioKanga has been developed by Dr Stuart Stephen, with contributions from other team member in CSIRO.
