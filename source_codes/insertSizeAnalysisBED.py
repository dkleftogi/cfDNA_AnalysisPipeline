#!/usr/local/bin/python

'''
				Analysis of coverage and fragment lenght for bed file with multiple regions
                            Parallel implementation using multi-threading

BEGIN COPYRIGHT NOTICE

    cfDNApipeline code -- (c) 2019 Dimitrios Kleftogiannis -- GIS -- A*STAR
    
    Copyright 2019 Genome Institute of Singapore (GIS) and Agency for Science, Technology and Research (A*STAR).

    This Program is free software licensed under the MIT License.

    You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENSE.md".

    This Program is distributed as a service to the research community and is experimental in nature and may have hazardous properties. 
    
    The Program is distributed WITHOUT ANY WARRANTY, express or implied. In particular all warranties as to SATISFACTORY QUALITY or FITNESS FOR A PARTICULAR PURPOSE are excluded. 

    Published reports of research using this code (or a modified version) should cite the relevant article of this tool.

    Comments and bug reports are welcome.
       
    Email to dimitrios.kleftogiannis@kaust.edu.sa

    I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
    You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

    
END COPYRIGHT NOTICE

UTILITY
  This program takes as input a bam file and a bed file with many intervals and computes
  per interval,the insert size distribution (TLEN field in the sam format)


INPUT ARGUMENTS
    
    1. bam file              : a bam file from the tumour or plasma of interest
    
    2. bed file              : a bed file with the genomic regions of interest

    3. output DIR             : the DIR (full path) to store the results

                                  
RUNNING
	
	An execution example is as follows:

    python insertSizeAnalysisBED.py bamFile=test.bam bedFile=intervals.bed outDIR=XXXXX

    To obtain the toy data we used for testing please contact Dimitrios

'''

import sys
import os
import pysam
import re
from collections import defaultdict
from itertools import groupby
import datetime
import time
import threading

#prints information about program's execution
def printUsage():
    print('To run this program please type the following:')
    print('\tpython insertSizeAnalysisBED bamFile=file.bam bedFile=file.bed outDIR=/your/output/dir\n')
    print('Where:\n') 
    print('\tfile.bam is a bam file. The bam file needs to be indexed with samtools.\n')
    print('\tfile.bed is a standard bed file with minimum three columns chrom start stop\n')
    print('\t/your/output/dir is the full path of the location to store the results\n')

    print('Please give the arguments in the indicated order similar to the provided example!\n') 

#global variables
exitFlag = 0
#define minimum mapping quality;
quality='10'

#class that processes one thread with input of one chromosome
class myThread (threading.Thread):
    def __init__(self, threadID, name, chrom,bamFilePrefix, bamFile,SAM_FILES):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.chrom = chrom
        self.fileName = bamFilePrefix
        self.inputBam = bamFile
        self.myDir = SAM_FILES
    def run(self):
        myStr="\t\tStarting " + self.name +" to process "+ self.chrom
        #print('%s\n'%myStr) 
        splitChrom(self.name,self.chrom,self.fileName, self.inputBam,self.myDir)

#run samtools for one chromosome
def splitChrom(threadName, chrom, fileName, inputBam, myDir):
   if exitFlag:
      threadName.exit()
   #print "%s: %s" % (threadName, time.ctime(time.time()))
   command='samtools view -b '+inputBam+' '+chrom+' > '+myDir+'/'+fileName+'_'+chrom+'.bam && samtools index '+myDir+'/'+fileName+'_'+chrom+'.bam'
   os.system(command)
   #print(command)
   #print('\t\tFinished:%s\n'%(threadName))

#store the bed file
def storeBedFile(bedFile):
    #simple dict to store the values
    aDict={}
    #check if file exists
    if os.path.exists(bedFile):
        #open the file and read it line by line
        fileIN=open(bedFile,'r')
        #read the first line and store it to the dictionary
        for eachLine in fileIN:
            line = eachLine.rstrip('\r')
            line = line.rstrip('\n')
            tmp=line.split("\t")
            chrom=tmp[0]
            startPos=tmp[1]
            endPos=tmp[2]
            key=chrom+'\t'+startPos+'\t'+endPos
            #save the record
            aDict[key]=chrom
        fileIN.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function storeBedFile: The bed file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()
    return aDict

#save the chromosomes
def storeChromosomes(bedFile):
    #save the positions
    aList=[]
    #check if bed file exists
    if os.path.exists(bedFile):
        #the file exists
        #open the file and read it line by line
        fileIN=open(bedFile,'r')
        #read the first line and store it to the dictionary
        for eachLine in fileIN:
            line = eachLine.rstrip('\n')
            tmp=line.split("\t")
            chrom=tmp[0]
            startPos=tmp[1]
            key=chrom
            aList.append(key)
            #print repr(key)
        fileIN.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function storeChromosomes: The bed file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()
    #remove duplicates from the list
    return list(set(aList))

#this function generates the report
def generateReport(bedList,bamFilePrefix,RESULTS,outDIR):

    #write the file
    reportFile=outDIR+'/'+bamFilePrefix+'_Distreport.txt'
    OutReportFile=open(reportFile,'w')
    #report a message
    count=1
    for idx in bedList:

        myChrom=bedList[idx]
        myInterval=idx
        mtInterval=myInterval.split('\n')
        #we do this just to index the output
        tmp=myInterval.split('\t')
        a=tmp[0]
        b=tmp[1]
        c=tmp[2]
        KEY=a+':'+b+'-'+c
        #generate the input sam file for the target region
        samFileName=RESULTS+'/'+bamFilePrefix+'_Target.sam'
        myBam=RESULTS+'/'+bamFilePrefix+'_'+myChrom+'.bam'
        generateSamFile(myBam,samFileName,myInterval,RESULTS)
        #read the target sam and generate the output
        parseSamFile(samFileName,OutReportFile,KEY)
        #uncomment the following code to monitor progress
        #if count%100==0:
         # print('Processed %d intervals\n'%(count))
        #count=count+1

    #close the output file
    OutReportFile.close()

#this function takes as input the bamFile, and the region of interest and produces a sam file with all reads spanning the positions
def generateSamFile(bamFile,samFileName,myInterval,RESULTS):
#first check if the bam file exists
    if os.path.exists(bamFile):
        #the bam file exists
        
        #write the interval in a bed file
        tmpBedFile=RESULTS+'/tmp.bed'
        tmpOut=open(tmpBedFile,'w')
        tmpOut.write('%s'%(myInterval))
        tmpOut.close()
        #prepare the output of SAM file
        outSAM=open(samFileName,'w')
        for eachLine in pysam.view('-q',quality,'-L',tmpBedFile,bamFile):
          outSAM.write(eachLine)
        outSAM.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function: generateSamFile. The bam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()


#this function parses the sam file and computes the number of reads as well as the number of 
def parseSamFile(samFileName,OutReportFile,KEY):
  if os.path.exists(samFileName):

    InSamFile=open(samFileName,'r')
    countTotalReads=0
    countReadPairs=0
    countNotPaired=0
    countNotPairedSame=0
    countNotPairedOther=0
    qHash=defaultdict(list)
    pairedLen=[]
    nonPairedLen=[]
    for eachLine in InSamFile:
      countTotalReads=countTotalReads+1
      #store the reads by QNAME to identify the paired reads
      line = eachLine.rstrip('\n')
      tmp=line.split("\t")
      #read the sam fields neeeded
      QNAME=tmp[0]
      FLAG=tmp[1]
      RNAME=tmp[2]
      POS=tmp[3]
      MAPQ=tmp[4]
      CIGAR=tmp[5]
      RNEXT=tmp[6]
      PNEXT=tmp[7]
      TLEN=tmp[8]
      #store records to the qHash
      #consider all reads irespective of the mapping quality
      qnameKey=QNAME
      qHash[qnameKey].append(eachLine)
    #parse the reads by Qname
    for key in qHash:
      pairFound=qHash[key]
      for idx in pairFound:
        #find the ones with no mate
        if len(pairFound)==1:
        #find the ones with read mate
          tmp=idx.split("\t")
          QNAME=tmp[0]
          FLAG=tmp[1]
          RNAME=tmp[2]
          POS=tmp[3]
          MAPQ=tmp[4]
          CIGAR=tmp[5]
          RNEXT=tmp[6]
          PNEXT=tmp[7]
          TLEN=tmp[8]
          myLen=abs(int(TLEN))
          #print(QNAME)
          countNotPaired=countNotPaired+1
          if RNEXT=='=':
            countNotPairedSame=countNotPairedSame+1
            nonPairedLen.append(myLen)
          else:
            #this is a translocation
            countNotPairedOther=countNotPairedOther+1
        elif len(pairFound)==2:
        #check for errors
          tmp=idx.split("\t")
          QNAME=tmp[0]
          FLAG=tmp[1]
          RNAME=tmp[2]
          POS=tmp[3]
          MAPQ=tmp[4]
          CIGAR=tmp[5]
          RNEXT=tmp[6]
          PNEXT=tmp[7]
          TLEN=tmp[8]
          myLen=int(TLEN)
          if myLen>0:
            pairedLen.append(myLen)
          countReadPairs=countReadPairs+1
        else:
          #used for debbuging
          #print('Malakia paizei:%s'%idx)
          myskip=0
          #used for debugging
    #print('Reads:%d\tReadsPaired:%d\tReadsNotPaired:%d\tReadsSame:%s\tReadsNotSame:%d\n'%(countTotalReads,countReadPairs,countNotPaired,countNotPairedSame,countNotPairedOther))
    #print('Paired len:%s\n\n\n'%pairedLen)
    #print('Non paired len:%s\n\n\n'%nonPairedLen)
    
    #OutReportFile.write('%s'%(KEY)),
    k=0
    for idx in pairedLen:
      if k==0:
        #OutReportFile.write('\t%d'%idx),
        OutReportFile.write('%d\n'%idx)
        k=k+1
      else:
        OutReportFile.write('%d\n'%idx)
        #OutReportFile.write('\t%d'%idx),
        k=k+1
    k=0
    for idx in nonPairedLen:
      if k==0:
        OutReportFile.write('%d\n'%idx)
        #OutReportFile.write('\t%d'%idx),
        k=k+1
      else:
        OutReportFile.write('%d\n'%idx)
        #OutReportFile.write('\t%d'%idx),
        k=k+1
    #OutReportFile.write('\n')
    del qHash
    del pairedLen
    del nonPairedLen
  else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function parseSamFile: The input sam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=4:
        printUsage()
    else:
        #parse the first input arguments 
        #if the user does not input the correct argument name it gets an error and the program stops
        bamFile=sys.argv[1].split('bamFile=')
        bamFile=bamFile[1]
        #parse the second argument
        bedFile=sys.argv[2].split('bedFile=')
        bedFile=bedFile[1]

        #parse the third argument
        outDIR=sys.argv[3].split('outDIR=')
        outDIR=outDIR[1]

        #get the bam file name
        bamFilePrefix=bamFile[:-4]
        tmp=bamFilePrefix.split('/')
        bamFilePrefix=tmp[-1]

       
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Function storeBedFile: store the genomic regions of interest'%(st))
        bedList=storeBedFile(bedFile)

        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Function storeChromosomes: store the chromosomes of interest'%(st))
        chromList=storeChromosomes(bedFile)
       
        #generate folders
        RESULTS=outDIR+'/'+bamFilePrefix
        command='mkdir -p '+RESULTS
        os.system(command)

        #split the original bam into one bam per chromosome
        #this is a multi-threading implementation that gives one thread per chromosome
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Spliting the original BAM per chromosome: Releasing threads...'%(st))
        #release the threads
        threadLock = threading.Lock()
        
        #wave 1 for load balancing
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr1' or chrom=='chrY' or chrom=='chr22' or chrom=='chr21' or chrom=='1' or chrom=='Y' or chrom=='22' or chrom=='21':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 1
        for t in threads:
            t.join()
        del threads

        #wave 2
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr2' or  chrom=='chr19' or chrom=='chr20' or chrom=='2' or chrom=='19' or chrom=='20':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 2
        for t in threads:
            t.join()
        del threads

        #wave 3
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr3' or  chrom=='chr18' or chrom=='chr17' or chrom=='3' or chrom=='18' or chrom=='17':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 3
        for t in threads:
            t.join()
        del threads

        #wave 4
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr4' or  chrom=='chr16' or chrom=='chr15' or chrom=='4' or chrom=='16' or chrom=='15':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 4
        for t in threads:
            t.join()
        del threads

        #wave 5
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr5' or  chrom=='chr14' or chrom=='chr13' or chrom=='5' or  chrom=='14' or chrom=='13':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 5
        for t in threads:
            t.join()
        del threads

        #wave 6
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr6' or  chrom=='chr12' or chrom=='chr11' or chrom=='6' or chrom=='12' or chrom=='11':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 6
        for t in threads:
            t.join()
        del threads

        #wave 7
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr7' or  chrom=='chr10' or chrom=='chr9' or chrom=='7' or chrom=='10' or chrom=='9':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 7
        for t in threads:
            t.join()
        del threads

        #wave 8
        threads = []
        i=1
        for chrom in chromList:

            if chrom=='chr8' or  chrom=='chrX' or chrom=='8' or  chrom=='X':
                myStr='Thread-'+str(i)
                thread = myThread(i, myStr, chrom,bamFilePrefix,bamFile,RESULTS)
                thread.start()
                threads.append(thread)
                i=i+1
        #collect wave 8
        for t in threads:
            t.join()
        del threads

        
        #collect them
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Spliting the original BAM per chromosome: All threads collected...'%(st))

        #generate the complete sam file
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Function generateReport: parse the bed file and generate the distribution of lenght'%(st))
        generateReport(bedList,bamFilePrefix,RESULTS,outDIR)

        #clean interm files
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Remove intermediate results to save some space....'%(st))
        command='rm '+RESULTS+'/*.bam && rm '+RESULTS+'/*.bai && rm '+RESULTS+'/*.sam && rm '+RESULTS+'/*.bed'
        os.system(command)
        
        print('************************************************************************************************************************************\n')

#start
if __name__=='__main__':
    myMain()

