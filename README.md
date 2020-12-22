# rEDIT

```
SYNOPSIS
        ./redit -q <path> -t <path> -l <length> -m (g|sg) -a (MEM|MUM) [--naive] [--load <prefix1>]
                [--save <prefix2>]

OPTIONS
        <path>      query sequences in fasta or fastq format
        <path>      target sequence in fasta format
        <length>    minimum match length (e.g., 20)
        g|sg        distance function (e.g., global or semi-global)
        MEM|MUM     filter maximal exact matches, 'MEM' consider all, 'MUM' must be single-copy
        --naive     use slow 2d dynamic programming algorithm for correctness check
        <prefix1>   load suffix array from files starting with specified prefix
        <prefix2>   save suffix array to files starting with specified prefix
```
