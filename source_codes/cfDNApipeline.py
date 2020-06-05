#!/usr/local/bin/python

'''
cfDNA analysis pipeline main script

BEGIN COPYRIGHT NOTICE

    cfDNApipeline code -- (c) 2019 Dimitrios Kleftogiannis -- GIS -- A*STAR
    
    Copyright 2019 Genome Institute of Singapore (GIS) and Agency for Science, Technology and Research (A*STAR).

    This Program is free software licensed under the MIT License.

    You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENSE.md".

    This Program is distributed as a service to the research community and is experimental in nature and may have hazardous properties. 
    
    The Program is distributed WITHOUT ANY WARRANTY, express or implied. In particular all warranties as to SATISFACTORY QUALITY or FITNESS FOR A PARTICULAR PURPOSE are excluded. 


    Published reports of research using this code (or a modified version) should cite the relevant article of this tool.

    Comments and bug reports are welcome.
       
    Email to dimitrios_kleftogiannis@gis.a-star.edu.sg

    I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
    You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY
  This program executes the cfDNA analysis workflow.


INPUT ARGUMENTS
    1. The directory you passed to myEnvConfig.sh script where the Conda env is installed

    2. An un-deDuped BAM file with UMI tags (absolute path)

    3. The directory to store the results (absolute path)

    4. A reference genome in fasta format (indexed)

    5. A bed file with the panel design (absolute path)

    6. The minimum variant allele frequency used for SNV calling (default is 0.005)

    7. The string tag incorporated to the BAM file, used to group reads from the same family (e.g., RX or ZU)

    8. The directory you downloaded the VEP Cache file (absolute path)

   
DEPENDENCIES
    To resolve the depending programs make sure you run the provided configuration program myEnvConfig.sh and you activated the conda env  
    

RUNNING
    
    An execution example is as follows:

    python cfDNApipeline.py workingDir=/home/centos fileName=/home/centos/data/example.bam resultsDir=/home/centos/Results referenceGenome=/home/centos/References/hg38.fa bedFile=/home/centos/bedFiles/panel_77.bed minVAF=0.005 tagUMI=ZU vepDir=/home/centos/vepCache


    To obtain toy data contact Dimitrios

'''


import sys
import os
import re
from itertools import groupby
from collections import defaultdict, Counter
import datetime
from datetime import date
import time
import threading


#prints information about program's execution
def printUsage():
    print('To run cfDNA analysis please follow the example below:\n')
    print('python cfDNApipeline.py workingDir=XX fileName=XX resultsDir=XX referenceGenome=XX bedFile=XX minVAF=XX tagUMI=XX vepDir=XX\n')
    print('Where:\n') 
    print('workingDir      is the absolute directory where the conda env is installed\n')
    print('fileName        is the absolute path of the bam file\n')
    print('resultsDir      is the absolute path to a directory for storage of results \n')
    print('referenceGenome is the reference genome in fasta format (absolute path)\n')
    print('bedFile         is the bed file with panel design (absolute path)\n')
    print('minVAF          is the minimum variant allele frequency used for variant calling (default 0.005)\n')
    print('tagUMI          is tag string used to group reads from the same family\n')
    print('vepDir          is the absolute directory where VEP Cache is downloaded and indexed\n')
    print('\n\nExample:\n')
    print('python cfDNApipeline.py workingDir=/your/working/dir fileName=/path/to/your/example.bam resultsDir=/your/results/dir referenceGenome=/path/to/your/reference/genome.fa bedFile=/path/to/your/panel.bed minVAF=0.005 tagUMI=ZU vepDir=/path/to/your/VEP/cache')
    print('\n\nPlease give the arguments in the indicated order similar to the example provided!\n') 
    print('\n\nAnd remember to activate the conda env\n') 


#find the bam file prefix and the actual path and save the results
def storeFile(myFile):
   #simple dict to store the values
    aDict={}
    #check if file exists
    if os.path.exists(myFile):
        
        line = myFile.rstrip('\n')
        tmp=line.split("/")
        sampleName=tmp[-1]
        sampleName=sampleName[0:-4]
        #print(filename)
        restFile='/'.join(tmp[0:-1])
        #print('%s with %s\n'%(sampleName,restFile))

        aDict[sampleName]=restFile
        
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function storeFile: The file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()
    return aDict

#this function takes the bam file with UMIs and removes the duplicates using  UMIs
def processSample(fileHash,workingDir,resultsDir,referenceGenome,scriptsFolder,tagUMI):

    #workingDir is the absolute path to your conda installation given from argument
    GATK_merge='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/picard-2.20.3-0/picard.jar MergeSamFiles'
    GATK_deDup='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/picard-2.20.3-0/picard.jar MarkDuplicates'
    GATK_sort='java  -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/picard-2.20.3-0/picard.jar SortSam'
    FGBIO_mateInfo='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/fgbio/fgbio.jar --tmp-dir=/tmp SetMateInformation'
    FGBIO_groupUMI='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/fgbio/fgbio.jar --tmp-dir=/tmp  GroupReadsByUmi'
    FGBIO_generateConsensus='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/fgbio/fgbio.jar --tmp-dir=/tmp  CallMolecularConsensusReads'
    FGBIO_filterConsensus='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/fgbio/fgbio.jar  -tmp-dir=/tmp  FilterConsensusReads'
    GATK_samToFastq='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/picard-2.20.3-0/picard.jar SamToFastq'
    GATK_mergeAlignment='java -Xmx10g -XX:-UseGCOverheadLimit -d64 -jar '+workingDir+'/miniconda2/share/picard-2.20.3-0/picard.jar MergeBamAlignment'
    BWA_MEM='bwa mem -t 8 -v 2 -R '


    for myArg in fileHash:

        outScriptFile=scriptsFolder+'/'+myArg+'_PrePro_commands.txt'
        outScript=open(outScriptFile,'w')

        projectName='cfDNApipelinePart1'
        mergedDIR=fileHash[myArg]

        #command 1

        #the last part of the command might change depending on the conda env --> e.g., in Amazon EC2 clusters provided by RONIN TMP_DIR=/shared is a dir with a lot of space to store interm results
        mainCommand=GATK_sort+' I='+mergedDIR+'/'+myArg+'.bam O='+resultsDir+'/'+myArg+'.mergedQNAMEsorted.bam SO=queryname VERBOSITY=INFO TMP_DIR=/tmp'
        outScript.write(mainCommand)
        outScript.write("\n\n")


        #command 2
        mainCommand=FGBIO_mateInfo+' -i='+resultsDir+'/'+myArg+'.mergedQNAMEsorted.bam -o='+resultsDir+'/'+myArg+'.mergedQNAMEsortedFixed.bam'
        outScript.write(mainCommand)
        outScript.write("\n\n")


        #command 3
        #the command take parameters --edits=1 --min-map-q=20 --strategy=adjacency  by default, this might by the user if needed (hardcoded)
        mainCommand=FGBIO_groupUMI+' --input='+resultsDir+'/'+myArg+'.mergedQNAMEsortedFixed.bam --output='+resultsDir+'/'+myArg+'.grouped.bam --edits=1 --min-map-q=20 --raw-tag='+tagUMI+' --strategy=adjacency --family-size-histogram='+resultsDir+'/'+myArg+'.groupedHist.txt'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #command 4
        mainCommand=GATK_sort+' I='+resultsDir+'/'+myArg+'.grouped.bam O='+resultsDir+'/'+myArg+'.groupedSorted.bam SO=coordinate TMP_DIR=/tmp'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #command 5
        mainCommand='samtools index '+resultsDir+'/'+myArg+'.groupedSorted.bam'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #command 6
        #the following parameters --error-rate-post-umi=30 --min-reads=2 --tag=MI are by default, this might by the user if needed (hardcoded)
        mainCommand=FGBIO_generateConsensus+' --input='+resultsDir+'/'+myArg+'.grouped.bam --output='+resultsDir+'/'+myArg+'.consensusUnMapped.bam --error-rate-post-umi=30 --min-reads=2 --tag=MI'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #command 7
        mainCommand=GATK_samToFastq+' I='+resultsDir+'/'+myArg+'.consensusUnMapped.bam F='+resultsDir+'/'+myArg+'.R1.fastq F2='+resultsDir+'/'+myArg+'.R2.fastq VALIDATION_STRINGENCY=SILENT TMP_DIR=/tmp'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #command 8
        rgTAG='\'@RG\\tID:FgBio.'+myArg+''+'\\tPL:ILLUMINA\\tLB:ALL\\tPU:NA\\tSM:'+projectName+'\\tCN:NA\''
        mainCommand=BWA_MEM+rgTAG+' -M '+referenceGenome+' '+resultsDir+'/'+myArg+'.R1.fastq '+resultsDir+'/'+myArg+'.R2.fastq | samtools view -hb - > '+resultsDir+'/'+myArg+'.UnSorted.FgbioDeDup.bam'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #command 9
        mainCommand=GATK_sort+' I='+resultsDir+'/'+myArg+'.UnSorted.FgbioDeDup.bam O='+resultsDir+'/'+myArg+'.FgbioDeDup.bam SO=coordinate TMP_DIR=/shared'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #command 10
        mainCommand='samtools index '+resultsDir+'/'+myArg+'.FgbioDeDup.bam'
        outScript.write(mainCommand)
        outScript.write("\n\n")

        outScript.close()

        #submit it
        command='sh '+outScriptFile
        os.system(command)

#this function runs variant screening and annotation
def variantScreening(fileHash,workingDir,resultsDir,referenceGenome,bedFile,scriptsFolder,minVAF_float,vepDir):

    for myArg in fileHash:

        outScriptFile=scriptsFolder+'/'+myArg+'_VarScreening_commands.txt'
        outScript=open(outScriptFile,'w')

        mainCommand='vardict -G '+referenceGenome+' -f '+str(minVAF_float)+ ' -N '+myArg+' -b '+resultsDir+'/'+myArg+'.FgbioDeDup.bam -z -c 1 -S 2 -E 3 -g 4 -h '+bedFile+' | '+workingDir+'/miniconda2/share/vardict-2019.06.04-0/teststrandbias.R | '+workingDir+'/miniconda2/share/vardict-2019.06.04-0/var2vcf_valid.pl -N '+myArg+' -E -f '+str(minVAF_float)+' > '+resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.vcf '
        outScript.write(mainCommand)
        outScript.write("\n\n")

        #to make it more complete we also annotate variants 
        #check here https://m.ensembl.org/info/docs/tools/vep/script/vep_cache.html for more info if needed
        mainCommand='vep -i '+resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.vcf -o '+resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.vcf --species homo_sapiens --cache --dir '+vepDir +' --canonical --check_existing --force_overwrite --vcf --buffer_size 50'
        outScript.write(mainCommand)
        outScript.write("\n")

        #the second script is ready and we just need to submit
        outScript.close()

        #submit it
        command='sh '+outScriptFile
        os.system(command)


#this function parses the annotated VCF and selects only the MODERATE/HIGH calls that fullfill specific criteria
def filterVCF(fileHash,resultsDir):

    for myArg in fileHash:
        myfile=resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.vcf'
        if os.path.exists(myfile):

            InVcfFile=open(myfile,'r')

            outFileName=resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.filtered.txt'
            outFile=open(outFileName,'w')

            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            outFile.write('#Original input VCF file: %s.FgbioDeDup.VarDict.VEP.vcf -- Report produced: %s\n'%(myArg,st))
            outFile.write('#CHROM\tPOS\tREF\tALT\tGENE\tvariantType\tCOV\tSupporintReads\tVAF\tRefBias\tVarBias\tMeanReadPos\tMeanReadQual\tMeanMappingQual\tSignalToNoise\tHighQualVarReads\tHighQualCov\tConsequence\tImpact\tBiotype\tProteinPos\tAminoAcidChange\tExistingVar\tPopulationVAF\tPredictedClinicalSignif\n')

            for eachLine in InVcfFile:

                #skipp the header
                if eachLine[0]!='#':

                    line=eachLine.rstrip('\n')
                    tmp=line.split('\t')
                    #parse the fields we need
                    CHROM = tmp[0]
                    POS = tmp[1]
                    ID = tmp[2]
                    REF = tmp[3]
                    ALT = tmp[4]
                    QUAL = tmp[5]
                    FILTER = tmp[6]

                    if FILTER=='PASS':
                        #now continue and get more specific fields about the variant of interest
                        FORMAT=tmp[7]
                        tmp=FORMAT.split(';')
                        #the last at the end is the VEP annotation that we need to parse separately using comma delimiter
                        annotationFields = tmp[-1]

                        sample = tmp[0].split('SAMPLE=')
                        sample = sample[1]

                        variantType = tmp[1].split('TYPE=')
                        variantType = variantType[1]

                        dp = tmp[2].split('DP=')
                        dp = int(dp[1])

                        vd = tmp[3].split('VD=')
                        vd = int(vd[1])

                        vaf = tmp[4].split('AF=')
                        vaf = float(vaf[1])

                        #not sure if we use this info
                        bias = tmp[5].split('BIAS=')
                        bias = bias[1]

                        refbias = tmp[6].split('REFBIAS=')
                        refbias = refbias[1]

                        varbias = tmp[7].split('VARBIAS=')
                        varbias = varbias[1]

                        #mean position in reads
                        pmean = tmp[8].split('PMEAN=')
                        pmean = float(pmean[1])
                        #mean quality in reads
                        qual = tmp[10].split('QUAL=')
                        qual = float(qual[1])
                        #this is the Fisher's test p-value from VarDict
                        sbf = tmp[12].split('SBF=')
                        sbf = float(sbf[1])

                        #mean mapping quality
                        mq = tmp[14].split('MQ=')
                        mq = float(mq[1])

                        #the higher the number the better the call
                        sn = tmp[15].split('SN=')
                        sn = float(sn[1])

                        #consider high coverage variant reads
                        hicnt = tmp[22].split('HICNT=')
                        hicnt = int(hicnt[1])

                        #consider high coverage variables
                        hicov = tmp[23].split('HICOV=')
                        hicov = int(hicov[1])

                        tmp = annotationFields.split(',')
                        for idx in tmp:
                            #find the canonical transcript
                            if 'YES' in idx:

                                canonicalTrans = idx
                                vepTmp = canonicalTrans.split('|')
                                Consequence = vepTmp[1]
                                IMPACT = vepTmp[2]
                                GENE = vepTmp[3]
                                BIOTYPE = vepTmp[7]
                                Protein_position = vepTmp[14]
                                Amino_acids = vepTmp[15]
                                Existing_variation = vepTmp[17]
                                AF = vepTmp[24]
                                CLIN_SIG = vepTmp[25]
                                #here we continue otherwise we skip because we are not interested in synonymous variants...
                                if IMPACT=='MODERATE' or IMPACT=='HIGH':
                                    #write the output except for Complex type of mutations returned by VarDict
                                    if 'Complex' not in variantType:
                                        
                                        #apply some filtering based on reads
                                        a=varbias.split(':')
                                        FW=int(a[0])
                                        BW=int(a[1])
                                        #here we consider only high quality reads, just to refine our results
                                        if hicov>=100 and hicnt>=3 and FW>1 and BW>1 and pmean>15 and sn>20:
                                            outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(CHROM,POS,REF,ALT,GENE,variantType,dp,vd,vaf,refbias,varbias,pmean,qual,mq,sn,hicnt,hicov,Consequence,IMPACT,BIOTYPE,Protein_position,Amino_acids,Existing_variation,AF,CLIN_SIG))
            InVcfFile.close()
            outFile.close()

        else:
            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('\n[%s] ERROR from function filterVCF: The input %s file does not exist!\n'%(st,myfile))
            print('************************************************************************************************************************************\n')
            sys.exit()

        
#this function parses the filtered VCF and generates the input required for duplexCaller to run
def convertFilteredVCF(fileHash,resultsDir):

    for myArg in fileHash:
        myfile=resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.filtered.txt'
        if os.path.exists(myfile):

            outputFile=resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.filtered_modified.txt'
            tmpFile=resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.filtered_tmp.txt'

            InFile=open(myfile,'r')
            outFile=open(tmpFile,'w')

            for eachLine in InFile:
                if eachLine[0]=='#':
                    a=1                    
                else:
                    line = eachLine.rstrip('\n')
                    tmp=line.split("\t")
                    chrom=tmp[0]
                    pos=tmp[1]
                    refAllele=tmp[2]
                    altAllele=tmp[3]
                    myType=tmp[5]
                    
                    if myType=='SNV' or myType=='SNP':
                        outFile.write('%s\t%s\n'%(chrom,pos))
                    elif myType=='DEL' or myType=='Deletion':
                        c=len(refAllele)
                        for myCount in range(int(pos),int(pos)+c):
                            outFile.write('%s\t%s\n'%(chrom,myCount))
                    elif myType=='INS' or myType=='Insertion':
                        c=len(altAllele)
                        for myCount in range(int(pos),int(pos)+c):
                            outFile.write('%s\t%s\n'%(chrom,myCount))
                    else:
                        ts = time.time()
                        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
                        print('\n[%s] WARNING from function produceOutput: No supported mutation type in %s!\n'%(st,inputFile))
                        print('************************************************************************************************************************************\n')
                        #sys.exit()
            InFile.close()
            outFile.close()

            command='sort '+tmpFile+' | uniq >'+outputFile
            os.system(command)
            #remove the tmp file
            command='rm '+tmpFile
            os.system(command)

        else:

            #the file does not exist
            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('\n[%s] ERROR from function convertFilteredVCF: The input file %s does not exist!\n'%(st,inputFile))
            print('************************************************************************************************************************************\n')
            sys.exit()


def runDuplexCaller(fileHash,resultsDir,referenceGenome,scriptsFolder):

    for myArg in fileHash:

        positionFile=resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.filtered_modified.txt'
        if os.path.exists(positionFile):

            outScriptFile=scriptsFolder+'/'+myArg+'_duplexCallerCommand.txt'
            outScript=open(outScriptFile,'w')
            myBAM=resultsDir+'/'+myArg+'.FgbioDeDup.bam'
            Num=1
            outScript.write('python duplexCallerModified_AWS.py bamFile=%s positionFile=%s referenceGenome=%s outDIR=%s index=%d\n'%(myBAM,positionFile,referenceGenome,resultsDir,Num))
            outScript.close()

            #we submit the duplexCaller script
            command='sh '+outScriptFile
            os.system(command)

        #the file does not exist
        else:
            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('\n[%s] ERROR from function runDuplexCaller: The input position file %s does not exist!\n'%(st,positionFile))
            print('************************************************************************************************************************************\n')
            sys.exit()


def cleanVariantReport(fileHash,resultsDir,referenceGenome,scriptsFolder):

    for myArg in fileHash:

        originalFile=resultsDir+'/'+myArg+'.FgbioDeDup.VarDict.VEP.filtered.txt'
        if os.path.exists(originalFile):

            variantReport=resultsDir+'/'+myArg+'.FgbioDeDup_VariantReport.txt'
            if os.path.exists(variantReport):

                outScriptFile=scriptsFolder+'/'+myArg+'_CleanVariantReportCommand.txt'
                outScript=open(outScriptFile,'w')
                outScript.write('python filterVariantReportAdjusted.py inputFile=%s originalVCF=%s\n'%(variantReport,originalFile))
                outScript.close()
                
                #submit it
                command='sh '+outScriptFile
                os.system(command)

            else:
                ts = time.time()
                st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
                print('\n[%s] ERROR from function cleanVariantReport: The variant report file %s does not exist!\n'%(st,variantReport))
                print('************************************************************************************************************************************\n')
                sys.exit()


        else:
            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('[%s] ERROR from function cleanVariantReport: The original file %s does not exist!\n'%(st,originalFile))
            print('************************************************************************************************************************************\n')
            sys.exit()


#this function generates the distribution of fragment length
def fragmentLenAnalysis(fileHash,resultsDir,bedFile,scriptsFolder):

    for myArg in fileHash:

        inputFile=resultsDir+'/'+myArg+'.FgbioDeDup.bam'
        if os.path.exists(inputFile):

            outScriptFile=scriptsFolder+'/'+myArg+'_fragLenCommand.txt'
            outScript=open(outScriptFile,'w')
            outScript.write('python insertSizeAnalysisBED.py bamFile=%s bedFile=%s outDIR=%s\n'%(inputFile,bedFile,resultsDir))
            outScript.close()

            #submit it
            command='sh '+outScriptFile
            os.system(command)

        else:

            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('[%s] ERROR from function fragmentLenAnalysis: The input bam file %s does not exist!\n'%(st,inputFile))
            print('************************************************************************************************************************************\n')
            sys.exit()

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=9:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\t\tYour input arguments are not correct!\n')
        print('\t\t\t\t\t   Genome Institute of Singapore (GIS) -- A*STAR\n')
        print('\t\t\tCopyright 2019 GIS -  Dimitrios Kleftogiannis & Jun Xian Liew - dimitrios_kleftogiannis@gis.a-star.edu.sg\n')
        #if arguments are not correct print a help message
        printUsage()
    else:
        print('************************************************************************************************************************************\n')
        print('\t\t\t cfDNApipeline.py: Run the full pipeline for cfDNA data processing and analysis\n')
        print('\t\t\t\t  Genome Institute of Singapore (GIS) -- A*STAR\n')
        print('\t\tCopyright 2019 GIS - Dimitrios Kleftogiannis & Jun Xian Liew - dimitrios_kleftogiannis@gis.a-star.edu.sg\n')        
        #parse the input arguments 
        
        workingDir   =sys.argv[1].split('workingDir=')
        workingDir   =workingDir  [1]

        fileName   =sys.argv[2].split('fileName=')
        fileName   =fileName  [1]
        
        resultsDir =sys.argv[3].split('resultsDir=')
        resultsDir =resultsDir[1]

        referenceGenome =sys.argv[4].split('referenceGenome=')
        referenceGenome =referenceGenome[1]

        bedFile =sys.argv[5].split('bedFile=')
        bedFile =bedFile[1]

        
        minVAF =sys.argv[6].split('minVAF=')
        minVAF =minVAF[1]
        
        minVAF_float=float(minVAF)
        if minVAF_float<0 or minVAF_float>1.00:
            print('\nWARNING: User gave invalid minVAF argument. Execution continues with default minVAF=0.005\n')
            minVAF_float=0.005

        tagUMI   =sys.argv[7].split('tagUMI=')
        tagUMI   =tagUMI  [1]

        vepDir   =sys.argv[8].split('vepDir=')
        vepDir   =vepDir  [1]

        #print the arguments given by user; is good for 'self' debugging
        print('Execution started with the following parameters:\n')
        print('1. workingDir             :         \t\t\t\t%s' % workingDir)
        print('2. fileName               :         \t\t\t\t%s' % fileName)
        print('3. resultsDir             :         \t\t\t\t%s' % resultsDir)
        print('4. referenceGenome        :         \t\t\t\t%s' % referenceGenome)
        print('5. bedFile                :         \t\t\t\t%s' % bedFile)
        print('6. minVAF                 :         \t\t\t\t%.4f' % minVAF_float)
        print('7. tagUMI                 :         \t\t\t\t%s' % tagUMI)
        print('8. vepDir                 :         \t\t\t\t%s' % vepDir)
       
        #generate a folder to store the scripts
        scriptsFolder=workingDir+'/cfDNApipeline_scripts'

        command='mkdir -p '+scriptsFolder #+' && mkdir -p '+LSF_logs+' && mkdir -p '+LSF_err
        os.system(command)
        #generate the folder for the results
        command='mkdir -p '+resultsDir
        os.system(command)

        #save the file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\n[%s] Function storeFile: store file name and path'%(st))
        fileHash=storeFile(fileName)
        
        #process the file and generate consensus using UMIs
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\n[%s] Function processSample: produce BAM file with consensus sequences'%(st))
        processSample(fileHash,workingDir,resultsDir,referenceGenome,scriptsFolder,tagUMI)

        #based on the consensus run variant screening
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\n[%s] Function variantScreening: perform variant screening and annotate variants'%(st))
        variantScreening(fileHash,workingDir,resultsDir,referenceGenome,bedFile,scriptsFolder,minVAF_float,vepDir)
        
        #filter the annotated VCF file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\n[%s] Function filterVCF: parse annotated VCF file and filter'%(st))
        filterVCF(fileHash,resultsDir)

        #prepare the data for duplexCaller
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\n[%s] Function convertInputVCF: parse filtered VCF file and convert position files'%(st))
        convertFilteredVCF(fileHash,resultsDir)

        #and run it
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\n[%s] Function runDuplexCaller: identify duplexes using duplexCaller'%(st))
        runDuplexCaller(fileHash,resultsDir,referenceGenome,scriptsFolder)

        #produce the filtered report
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function cleanVariantReport: produce final variant report'%(st))
        cleanVariantReport(fileHash,resultsDir,referenceGenome,scriptsFolder)

        #finally we run the fragment lenght analysis script insertSizeAnalysisBED.py
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function fragmentLenAnalysis: perform analysis of fragment lenght'%(st))
        fragmentLenAnalysis(fileHash,resultsDir,bedFile,scriptsFolder)
        

        print('************************************************************************************************************************************\n')

#this is where we start
if __name__=='__main__':
    myMain()
