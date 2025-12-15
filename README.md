# ChainX

ChainX is a tool that computes co-linear chaining costs between an input target and query sequences. It supports global and semi-global comparison modes. A unique aspect of ChainX is that it supports anchor overlaps and gap costs. The output chains are either optimal or close to optimal in practice.  

For a pair of sequences, computing chaining cost can be orders of magnitude faster than computing edit distance. Moreover, chaining cost and edit distance correlate well with each other. As a result, ChainX can serve as a faster alternative to estimating edit distance. More details about the algorithm are available in our [paper](#pub).

## Dependencies / External libraries
ChainX repository uses many third-party libraries. These are separately provided in [ext](ext) folder. 

- A C++ compiler with c++11 support, e.g., GNU g++ (version 5+)
- [essaMEM](https://doi.org/10.1093/bioinformatics/btt042)
- [libdivsufsort](https://github.com/y-256/libdivsufsort)
- [MUMmer](https://github.com/mummer4/mummer)
- [edlib](https://github.com/Martinsos/edlib)
- [clipp](https://github.com/muellan/clipp)
- [cxx-prettyprint](https://github.com/louisdx/cxx-prettyprint)
- [kseq](https://github.com/lh3/seqtk)

## Installation
```sh
git clone https://github.com/AT-CG/ChainX.git
cd ChainX
make
```

## Usage
```
SYNOPSIS
        ./chainX [-l <length>] [-a (MEM|MUM)] [--all2all] [--naive] -m (g|sg) -q <qpath> -t <tpath>

OPTIONS
        <length>    minimum anchor match length (default = 20)
        MEM|MUM     anchor type (default = MUM)
        --all2all   output all to all global distances among query sequences in phylip format
        --naive     use slow 2d dynamic programming algorithm to obtain exact cost
        g|sg        distance function (e.g., global or semi-global)
        <qpath>     query sequences in fasta or fastq format
        <tpath>     target sequence in fasta format
```

## Example
Test data can be accessed from [data](data) folder. Here is an example run.

```
$ ./chainX -m g -q data/time_global/mutated_80_perc.fasta -t data/time_global/Chromosome_2890043_3890042_0.fasta
INFO, chainx::parseandSave, target sequence file = data/time_global/Chromosome_2890043_3890042_0.fasta
INFO, chainx::parseandSave, query sequences file = data/time_global/mutated_80_perc.fasta
INFO, chainx::parseandSave, mode = g
INFO, chainx::parseandSave, anchor : minimim length = 20, type = MUM
INFO, chainx::main, read 1 queries, 999760 residues
INFO, chainx::main, read target, 1000000 residues
INFO, chainx::main, timer set
INFO, chainx::main, suffix array computed in 0.162169 seconds
INFO, chainx::main, timer reset
INFO, chainx::main, count of anchors (including dummy) = 2964, average length = 24.0651
INFO, chainx::main, query #0 (999760 residues), distance = 939635
INFO, chainx::main, distance computation finished (0.197675 seconds elapsed)
```

## <a name="pub"></a>Publications

- **Chirag Jain, Daniel Gibney and Sharma Thankachan**. "Algorithms for Colinear Chaining with Overlaps and Gap Costs". *Journal of Computational Biology (RECOMB 2022 special issue)*. [PDF](https://cds.iisc.ac.in/faculty/chirag/pubs/2022-jain-chainX-jcb.pdf)
