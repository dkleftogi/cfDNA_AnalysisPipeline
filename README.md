# Detection of genomic alterations using error-corrected cfDNA sequencing

Several recent studies have demonstrated the ability of cfDNA sequencing to provide early prognostication, better molecular profiling and monitoring of disease dynamics with many applications in genomic-driven oncology. 

We provide a bioinformatics pipeline that offers:
- cfDNA data preprocessing using UMIs
- ultra-sensitive detection of SNVs using VarDict and duplexCaller 
- annotation of variants using VEP and evaluation of fragmentation profiles of cfDNA

## Publication

Title: Detection of genomic alterations in breast cancer using circulating free DNA sequencing  

Journal: N/A

Published: pre-print available at bioRxiv

DOI: 

## Dependencies and System Requirements

Our bioinformatics pipeline is developed with Python version 2.7.10 and requires the following packages:

- bwa version 0.7.15
- ensembl-vep version 96.0
- fgbio version 0.8.1
- picard version 2.20.3
- pysam version 0.15.0.1
- R version 3.6.0
- samtools version 1.9
- vardict version 2019.06.04

To resolve all of the abovementioned dependenices we recommend installing these packages in a conda environment. You may refer to the following bash script which installs miniconda and all the required packages automatically.

```
myEnvConfig.sh
```

Please refer to Execution_examples.md for more information. We also provide step by step execution guidelines.

Our bioinformatics pipeline is developed on a Mac OS computer with Mojave version 10.14.5 and tested on Amazon EC2 m4.xlarge machines with Centos 7 operating system installed.

#### Current Release

19-Jul-2019 : Beta version 1 (installed and tested on Centos); does not clean intermediate results

## Contact

Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios_kleftogiannis@gis.a-star.edu.sg) OR Liew Jun Xian (liewjx@gis.a-star.edu.sg)

We are also interested to know about how you have used our source code, including any improvements that you have implemented.
 
You are free to modify, extend or distribute our source code, as long as our copyright notice remains unchanged and included in its entirety. 

## Licence

This project is licenced under the MIT Licence.

Copyright 2019 Genome Institute of Singapore (GIS) and Agency for Science, Technology and Research (A*STAR).

You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENCE".
