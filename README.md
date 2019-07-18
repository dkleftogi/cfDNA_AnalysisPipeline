# Detection of genomic alterations using error-corrected cfDNA sequencing

Several recent studies have demonstrated the ability of cfDNA sequencing to provide early prognostication, better molecular profiling and monitoring of disease dynamics with many applications in genomic-driven oncology. 

We provide a bioinformatics pipeline that offers ctDNA data preprocessing using UMIs, ultra-sensitive detection of SNVs using duplexCaller (available at https://github.com/dkleftogi/duplexFiltering) and evaluation of fragmentation profiles of cfDNA.


## Publication

Title: Detection of genomic alterations in breast cancer using circulating free DNA sequencing  

Journal: N/A

Published: pre-print available at bioRxiv

DOI: 

## Dependencies and System Requirements

The developed bioinformatics pipeline is written in Python 2 version 2.7.10.

The pipeline depends on the following packages:

```
fgbio version 0.8.1
```

```
bwa version 0.7.15
```

```
picard version 2.20.3
```

```
samtools version 1.9
```

```
pysam version 0.15.0.1
```

```
vardict version 2019.06.04
```

```
vardict version 2019.06.04
```

```
ensembl-vep vep 96.0
```

To resolve all the above mentioned dependenices we use Conda and we provide a bach script (myEnvConfig.sh ) that installs all required packages using Bioconda.


The pipeline has been developed in a Mac OS computer with Mojave version 10.14.5 and tested in AWS systems XXX and XXXX


#### Current Release

17-Jul-2019 : 


## Contact

Dr Dimitrios Kleftogiannis and Jun Xian Liew

Comments and bug reports are welcome, email to dimitrios_kleftogiannis@gis.a-star.edu.sg OR iewjx@gis.a-star.edu.sg

I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
You are free to modify, extend or distribute this code, as long as our copyright notice is included whole and unchanged. 

## Licence

Copyright 2019 -- Dimitrios Kleftogiannis -- Genome Institute of Singapore (GIS)

Agency of Science Research and Technology (A*STAR)
       			
You may not use this program except in compliance with the License. You may obtain a copy of the License at: https://opensource.org/licenses/ECL-2.0


