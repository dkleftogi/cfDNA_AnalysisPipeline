## FAQ

### Question 1: How do I install the pipeline ?

We provide a bash script that resolves all dependencies using Bioconda.

Run:

```
	sh myEnvConfig.sh myDIR

	where myDIR is the directory to download and generate the Conda env. 
```

Here is an example execution:

```
sh myEnvConfig.sh /home/centos 
```

PLEASE note that VEP requires a Cache to be downloaded and indexed. This step is usually slow and is not included in our configuration script. 

To install the Cache type the following commands once myEnvConfig.sh installation is COMPLETED.

```
mkdir /home/centos/vepCache

cd /home/centos/vepCache

curl ftp://ftp.ensembl.org/pub/release-96/variation/indexed_vep_cache/homo_sapiens_vep_96_GRCh38.tar.gz -o homo_sapiens_vep_96_GRCh38.tar.gz

tar xzf homo_sapiens_vep_96_GRCh38.tar.gz
```

More info about Ensembl VEP Caches can be found at https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html


### Question 2: How do I run the pipeline ?

We provide the following scripts:

```
cfDNApipeline.py

duplexCallerModified_AWS.py

filterVariantReportAdjusted.py

insertSizeAnalysisBED.py
```

that implement our analysis workflow. 

The main program is cfDNApipeline.py that take as input the following arguments:

```
	1. The directory you passed to myEnvConfig.sh script where the Conda env is installed

	2. An un-deDuped BAM file with UMI tags (absolute path)

	3. The directory to store the results (absolute path)

	4. A reference genome in fasta format (indexed)

	5. A bed file with the panel design (absolute path)

	6. The minimum variant allele frequency used for SNV calling (default is 0.005)

	7. The string tag incorporated to the BAM file, used to group reads from the same family (e.g., RX or ZU)

	8. The directory you downloaded the VEP Cache file (absolute path)
	
```

To make the example more interactive let us assume the following input arguments

```
	1. /home/centos 

	2. /home/centos/data/example.bam

	3. /home/centos/Results

	4. /home/centos/References/hg38.fa

	5. /home/centos/bedFiles/panel_77.bed

	6. 0.005

	7. ZU

	8. /home/centos/vepCache
	
```

In this scenario (first make sure the Conda env is activated) and then run:

```
source /home/centos/miniconda2/bin/activate

python cfDNApipeline.py workingDir=/home/centos fileName=/home/centos/data/example.bam resultsDir=/home/centos/Results referenceGenome=/home/centos/References/hg38.fa bedFile=/home/centos/bedFiles/panel_77.bed minVAF=0.005 tagUMI=ZU vepDir=/home/centos/vepCache

```

Please write the argyments in this order and be sure that the absolute paths are correct. We reccomend to use absolute paths, since this gives flexibility on cloud infrastructures. 


### Question 3: What the pipeline is doing and what are the outputs ?

If you execute the previous example correctly


