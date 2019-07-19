## FAQ

### Question 1: How do I resolve the dependencies ?

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

	8. The directory you downloaded the VEP Cache file (absolute path, see the previous comment about Caches).
	
```

To make the example more interactive let us assume the following exxample with input arguments:

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

In this scenario (always make sure the Conda env is activated, if not type source /home/centos/miniconda2/bin/activate) 

run:

```

python cfDNApipeline.py workingDir=/home/centos fileName=/home/centos/data/example.bam resultsDir=/home/centos/Results referenceGenome=/home/centos/References/hg38.fa bedFile=/home/centos/bedFiles/panel_77.bed minVAF=0.005 tagUMI=ZU vepDir=/home/centos/vepCache

```

Please follow the example, and write the arguments in this order. Always be sure that the absolute paths are correct. We reccomend to use absolute paths, since this gives flexibility on cloud infrastructures. 


### Question 3: What the pipeline is doing and what are the outputs ?

If you execute the previous example correctly, our pipeline processes file /home/centos/data/example.bam and performes UMI-aware deduplication using fgbio, picard and samtools. 

The output is a file named:
```
 /home/centos/Results/example.FgbioDeDup.bam
```
which contains only consensus reads after merging together reads with the same UMI tag.  

Please note that this step generates many intermediate results that are removed at the end of execution (TODO). 

Next, file /home/centos/Results/example.FgbioDeDup.bam is used as input to VarDict for variant screening with minimum VAF=0.005

The output of this step is a file named:
```
 /home/centos/Results/example.FgbioDeDup.Vardict.vcf
```

Then, the file /home/centos/Results/example.FgbioDeDup.Vardict.vcf is passed to VEP for annotation using the locally installed Cache (give the correct DIR as input).

The output of this step is a file named:
```
 /home/centos/Results/example.FgbioDeDup.Vardict.VEP.vcf
```

Next, we post-filter the annotated VCF file and select only MODERATE and HIGH impact variants that meet certain quality criteria as described in our manuscript. 

The output of this step is a file named:
```
 /home/centos/Results/example.FgbioDeDup.Vardict.VEP.filtered.txt
```

The post-filtered file is then passed to duplexCaller that identifies variants supported by in-silico duplexes. 

The original duplexCaller implementation can be found at https://github.com/dkleftogi/duplexFiltering

However here we provide a platform-independent version that over-simplifies execution.

The final output of variant calling is a file named:
```
 /home/centos/Results/example.FgbioDeDup_VariantsClean.txt
```

The last step of the pipeline concerns the evaluation profiles of cfDNA.

To do so we have developed a program named insertSizeAnalysisBED.py that evaluates fragmentation profiles of cfDNA and produces the following:
```
 /home/centos/Results/example.FgbioDeDup_Distreport.txt
```

### Question 4: Is it possible to run specific parts of the workflow and not the full thing ?

The short answer to this question is NO. We are planning to provide this functionality in the future.

However for the time being we reccomend to run alone only insertSizeAnalysisBED.py

This program takes as input the following argumenents:
```
	1. A BAM file with or without UMIs (absolute path)

	2. A bed file with the panel design (absolute path)

	3. The directory to store the results (absolute path)
	
```
Please note can produce the same results as part of the main cfDNApipeline.py program or using the following command:
```
python insertSizeAnalysisBED.py bamFile=/home/centos/data/TYL0002.bam bedFile=/home/centos/PanelDesigns/example_panel.bed outDIR=/home/centos/Results
```

### Question 5: I found a folder named cfDNApipeline_scripts, what does it contain ?

The folder is generated internally to store the execution commands. 

The commands are grouped based on the steps of the analysis as follows:

```
	1. example.FgbioDeDup_PrePro_commands.txt contains all commands for data pre-processing and deDup with UMIs
	2. example.FgbioDeDup_VarScreening_commands.txt contains all commands for variant screening and annotation 
	3. example.FgbioDeDup_duplexCallerCommand.txt contains the commands to execute duplexCaller
	4. example.FgbioDeDup_CleanVariantReportCommand.txt is the command the produces the final variant report
	5. example.FgbioDeDup_fragLenCommand.txt is the command to perform fragment lenght analysis

```
In case one of the steps fails, or you need to re-run some part you can simply re-submit the corresponding script. This saves time and effort.


### Question 6: I have problems to run the analysis, what should I do ?

Please read carefully the error messages and try to resolve any problem with dependencies. There might be also some problem with the paths, so be careful with your arguments.

Try to repeat the pipeline following the commands provided in the examples.

You can always contact Dimitrios or Jun Xian and ask for help =)


