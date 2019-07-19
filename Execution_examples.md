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

curl ftp://ftp.ensembl.org/pub/release-97/variation/indexed_vep_cache/homo_sapiens_vep_96_GRCh38.tar.gz -o homo_sapiens_vep_96_GRCh38.tar.gz

tar xzf homo_sapiens_vep_96_GRCh38.tar.gz
```

More info about Ensembl VEP Caches can be found at https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html







