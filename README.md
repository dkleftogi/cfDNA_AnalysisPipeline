# Detection of genomic alterations using error-corrected cfDNA sequencing

Several recent studies have demonstrated the ability of cfDNA sequencing to provide early prognostication, better molecular profiling and monitoring of disease dynamics with many applications in genomic-driven oncology. 

We provide a bioinformatics pipeline that offers ctDNA data preprocessing using UMIs, ultra-sensitive detection of SNVs using VarDict and duplexCaller, annotation of variants using VEP and evaluation of fragmentation profiles of cfDNA.


## Publication

Title: Detection of genomic alterations in breast cancer using circulating free DNA sequencing  

Journal: N/A

Published: pre-print available at bioRxiv

DOI: 

## Dependencies and System Requirements

The developed bioinformatics pipeline is written in Python 2 version 2.7.10.

The pipeline depends on the following packages:

fgbio version 0.8.1 ; bwa version 0.7.15 ; picard version 2.20.3 ; samtools version 1.9 ; pysam version 0.15.0.1 ; vardict version 2019.06.04 ; ensembl-vep version 96.0 ; R version 3.6.0

To resolve all the above mentioned dependenices we use Conda and we provide a bash script named 

```
myEnvConfig.sh
```

that installs automatically all required packages using Bioconda. 

Please see the Execution_examples.md for more information. We also provide step by step execution guideliness.


The pipeline has been developed in a Mac OS computer with Mojave version 10.14.5 and tested in AWS M4.XLARGE machines with Centos 7 and XXXX


#### Current Release

19-Jul-2019 : Beta version 1 (installed and tested on Centos); does not clean intermediate results


## Contact

Dimitrios Kleftogiannis and Liew Jun Xian

Comments and bug reports are welcome, email to dimitrios_kleftogiannis@gis.a-star.edu.sg OR liewjx@gis.a-star.edu.sg

I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
You are free to modify, extend or distribute this code, as long as our copyright notice is included whole and unchanged. 

## Licence

This project is licenced under the GNU GPLv3 Licence.

Copyright (C) 2019 -- Dimitrios Kleftogiannis -- Genome Institute of Singapore (GIS)

Agency of Science Research and Technology (A*STAR)
       			
You may not use this file except in compliance with the License. A copy of the licence is available, see file LICENCE.md 


