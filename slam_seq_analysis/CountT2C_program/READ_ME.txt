# slam-seq

A programme that processes alignments and detects T->C conversions
that is useful to process slam seq data very fast. It outputs text files.

Contains the library seqan in version 2.4.0 (for this library:
Copyright (c) 2006-2018, Knut Reinert, FU Berlin). See seqan/LICENCE
and https://www.seqan.de/). Does not work with newer versions of seqan.

To compile, you will need a c++ compiler compatible with c++14

Subdirectory CountT2C provides a programme that works with bulk RNA
seq data

Subdirectory EstPconv contains a script to estimate the T->C
conversion rate for a library

Subdirectory CountT2C_10X contains a programme similar to CountT2C
that handles 10X single cell sequencing data (preserves UMI and
Cell barcodes)

After compilation, to get help please use
./countT2C --help 





