kamix - Query k-mers matrices index as BGZF files.
======================================================

``kamix`` leverages the fantastic BGZF library in ``samtools`` to provide random access into
k-mer matrices files that have been compressed with ``bgzip``.

Kamix expect k-mer matrices compressed with bgzip. In this matrice, k-mer must be sorted in lexicographic.

Example :

```
tag	SRR2966453	SRR2966456	SRR2966471	SRR2966474	SRR2966454	SRR2966457	SRR2966472	SRR2966475	SRR2966455	SRR2966458	SRR2966473	SRR2966476
AAAAAATGTTTTGTAAGAAT	3	7	4	0	3	0	0	8	0	5	0	6
AAAAAATGTTTTGTAAGGAC	3	0	0	0	0	0	0	0	0	0	0	0
AAAAAATGTTTTGTAATTGA	5	5	5	7	3	8	4	5	6	0	3	9
AAAAAATGTTTTGTACAAAA	8	6	6	4	5	0	6	6	7	0	6	4
AAAAAATGTTTTGTAGAAAC	0	3	0	0	0	0	0	0	0	0	0	0
AAAAAATGTTTTGTAGAAAT	0	0	6	8	0	0	3	5	0	5	6	0
AAAAAATGTTTTGTAGACAT	6	4	0	6	3	0	0	0	3	0	0	0
```

## Examples

```
# 1. Index the k-mer matrice with bgzipped
bgzip counts-matrix.tsv

# 2. create a kamix index (big.vcf.gz.gbi)
kamix index counts-matrix.tsv.gz

# 3. query a k-mer." << endl;
kamix query counts-matrix.tsv.gz AAAAAAAAAAGGCTAAACAT

# extract the 100 random kmers.
kamix random counts-matrix.tsv.gz 100

# Is the file bgzipped?
kamix check counts-matrix.tsv.gz

# get total number of lines in the file (minus the header).
kamix size counts-matrix.tsv.gz
```


## Create a k-mer matrice with jellyfish and JoinCounts

**Counting and sorting 32-mer with DSK or Jellyfish**

* with DSK:

```
dsk -file sample.fastq.gz
dsk2ascii -file sample.h5 -out >(sort -k 1) > counts.tsv
```

* with Jellyfish:

```
jellyfish count -m 32 -s 10000 -o sample.jf <(zcat sample.fastq.gz)
jellyfish dump  -c sample.jf | sort -k 1 > counts.tsv
```

Consider using the sort command with `-S {resources.ram}G` and `--parallel {threads}` parameters to speed-up the sorting for large k-mer libraries.

**Join counts from multiple libraries with joinCounts**

Download and install [joinCounts](https://github.com/Transipedia/dekupl-joinCounts), to join counts files generated for each libraries.

```
joinCounts counts1.tsv counts2.tsv > counts-matrix.tsv
```

## Credit

kamix is derived from [grabix](https://github.com/arq5x/grabix) that uses BGZF library.
