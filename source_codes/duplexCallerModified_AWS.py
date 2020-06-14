#!/usr/local/bin/python

'''
duplexCaller version for identification of variants supported by reads with duplex configuration in deep-sequencing data

This code is originally published at https://github.com/dkleftogi/duplexFiltering and here we re-implemented it as part of a cfDNA analysis pipeline.

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
  This program takes as input a bam file and a bed-like file (two columns chrom TAB position) with genomic positions 
  of interest and applies specific filtering criteria for accepting or rejecting the variant call.


INPUT ARGUMENTS
    
    1. bam file                     : a bam file from the tumour or plasma of interest
    
    2. bed-like file                : a bed-like file with the positions of interest (chrom TAB position)
                              
    3. reference                    : the reference genome in fasta format

    4. output DIR                   : the DIR (full path) to store the results
    
    5. index                        : give the index of the chunk

RUNNING
	
	A simple execution example is as follows:

    python duplexCallerModified_AWS.py bamFile=test.chr17.bam positionFile=positions.bed referenceGenome=file.fa outDIR=Results index=1

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
    print('\tpython duplexCallerModified_AWS.py bamFile=file.bam positionFile=file.bed referenceGenome=file.fa outDIR=Results index=1\n')
    print('Where:\n') 
    print('\tfile.bam is a bam file.\n')
    print('\tfile.bed is a bed-like file with your input genomic positions. Please write chromosome and start only\n')
    print('\tfile.fa  is a reference genome in fasta format. Please use the same reference you used for the aligment of the bam file.\n')
    print('\tResults is a dir with full path to store the results\n')
    print('\tNum is the index of the chunk to process')
    print('\tThe reference genome and the input bam need to be indexed (produce files .bai, .fai and .dict) \n')
    print('Please give the arguments in the indicated order similar to the provided example!\n') 

exitFlag = 0

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
      
class myQnameThread (threading.Thread):
    def __init__(self, threadID, name, chrom,bamFilePrefix, bamFile,SAM_FILES,qHash):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.chrom = chrom
        self.fileName = bamFilePrefix
        self.inputBam = bamFile
        self.myDir = SAM_FILES
        self.qHash = qHash
    def run(self):
        myStr="\t\tStarting " + self.name +" to process "+ self.chrom
        #print('%s\n'%myStr) 
        parseQname(self.name,self.chrom,self.fileName, self.inputBam,self.myDir,self.qHash)

class myVariantThread (threading.Thread):
    def __init__(self, threadID, name, posKey,targetSamFileName,SAM_FILES):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.pos = posKey
        self.targetSamFileName = targetSamFileName
        self.myDir = SAM_FILES
    def run(self):
        myStr="\t\tStarting " + self.name +" to process position "+ self.pos
        #print('%s\n'%myStr) 
        splitFlatFile(self.name,self.pos,self.targetSamFileName,self.myDir)

#this is the function that parses the slice of the original QNAME and writes them 
def parseQname(threadName, chrom, fileName, inputBam, myDir,qHash):
    if exitFlag:
        threadName.exit()
    qnameFileName=myDir+'/'+fileName+'_'+chrom+'_QNAME_tmp.txt'
    fileOUT=open(qnameFileName,'w')
    for i in qHash:
        line = i.rstrip('\n')
        tmp=line.split("\t")
        QNAME=tmp[0]
        fileOUT.write('%s\n'%QNAME)
    fileOUT.close()
    command='sort '+qnameFileName+' | uniq >> '+myDir+'/'+fileName+'_'+chrom+'_QNAME.sam;'
    #print(command)
    os.system(command)

#this is the actual command
def splitChrom(threadName, chrom, fileName, inputBam, myDir):
   if exitFlag:
      threadName.exit()
   command='samtools view -b '+inputBam+' '+chrom+' > '+myDir+'/'+fileName+'_'+chrom+'.bam && samtools index '+myDir+'/'+fileName+'_'+chrom+'.bam'
   os.system(command)
   #print(command)
   #print('\t\tFinished:%s\n'%(threadName))

#save the genomic positions of interest; remember this is not the actula panel design but a modified version
def storePositionsFile(bedFile):
    #save the positions
    aDict={}
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
            key=chrom+'_'+startPos
            aDict[key]=int(startPos)
            #print repr(key)
        fileIN.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function storePositionsFile: The bed file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()
    return aDict

#this function takes as input the bamFile, the bedFile and produces a sam file with all reads spanning the positions
def generateTargetSamFile(posList,bamFile,targetSamFileName,SAM_FILES,index):
#first check if the bam file exists
    bamFilePrefix=bamFile[:-4]
    tmp=bamFilePrefix.split('/')
    bamFilePrefix=tmp[-1]
    if os.path.exists(bamFile):
        tmpSamFile=SAM_FILES+'/'+bamFilePrefix+'_tmpTarget'+index+'.sam'
        outSamFile=open(tmpSamFile,'w')
        tmpPosFile=SAM_FILES+'/posDict'+index+'.txt'
        outPosFile=open(tmpPosFile,'w')
        for key in posList:
            line = key.rstrip('\n')
            tmp=line.split('_')
            chrom=tmp[0]
            position=tmp[1]
            tmpBedFile=SAM_FILES+'/posTMP'+index+'.bed'
            tmpOut=open(tmpBedFile,'w')
            tmpOut.write('%s\t%s\n'%(chrom,position))
            tmpOut.close()
            tmpBamFile=SAM_FILES+'/'+bamFilePrefix+'_'+chrom+'.bam'
            #print('\t\tProcessing:%s\n'%(chrom))
            for eachLine in pysam.view('-F','256','-L',tmpBedFile,tmpBamFile):
                if eachLine=='\n':
                    outPosFile.write('%s\n'%(key))
                outSamFile.write(eachLine)

        outSamFile.close()
        outPosFile.close()
        #now generate the dictionary on disc...
        command='paste '+tmpPosFile+' '+tmpSamFile+' > '+ targetSamFileName
        os.system(command)

        command='rm '+tmpPosFile
        os.system(command)
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function generateTargetSamFile: The bam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()


#save the genomic positions of interest; this is not the panel design but a modified version
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

#parse the sam file and generate the flat file
def generateFlatFile(samFileName,flatFileName,refGenome):

    if os.path.exists(samFileName):
        InSamFile=open(samFileName,'r')
        #mark the 'bad cigars'
        countH=0
        countP=0
        if os.path.exists(refGenome):
            OutFile=open(flatFileName,'w')
            #Flag Variant is YES if the base in the position is similar to the reference, or NO otherwise. 
            OutFile.write('#INPUT_POS_ID\tQNAME\tCHROM\tREAD_MAPQ\tREAD_POS\tREF_POS\tCIGAR_FLAG\tREF\tCURRENT_BASE\tVARIANT\tQUAL\tRNAME:POS:RNEXT:PNEXT:TLEN:XI\n')
            #load the reference genome
            myFasta=pysam.FastaFile(refGenome)
            for eachLine in InSamFile:
                #here might be is a problem with the TAB delimiting in Python 3; check python or OS versions if there is any problem
                line = eachLine.rstrip('\n')
                tmp=line.split("\t")
                #read the sam fields we neeed
                posKey=tmp[0]
                QNAME=tmp[1]
                FLAG=tmp[2]
                RNAME=tmp[3]
                POS=tmp[4]
                MAPQ=tmp[5]
                CIGAR=tmp[6]
                RNEXT=tmp[7]
                PNEXT=tmp[8]
                TLEN=tmp[9]
                SEQ=tmp[10]
                QUAL=tmp[11]
                
                #skip the reads with low mapping quality
                MAPQ_int=int(MAPQ)
                if MAPQ_int>5:
                    #mask XI field
                    t='2'
                    extraField=RNAME+':'+POS+':'+RNEXT+':'+PNEXT+':'+TLEN+':'+t+':'+FLAG
                    numValues=re.findall('\d+',CIGAR)
                    flagValues=re.findall('[A-Za-z]',CIGAR)
                    TLEN_int=abs(int(TLEN))

                    prompt_pos=posKey.split('_')
                    prompt_pos=int(prompt_pos[1])
                    span_ref=0
                    #here we find all fragments with adapter read-through
                    for eachValue1,eachFlag1 in zip(numValues,flagValues):
                        if eachFlag1=='M':
                            if TLEN_int<int(eachValue1):
                                #find the real overlap of the reads
                                a=int(POS)
                                b=a+int(eachValue1)-1
                                c=int(PNEXT)
                                d=int(PNEXT)+TLEN_int-1
                                overlap_start=max(a,c)
                                overlap_end=min(b,d)
                                #the position is OK so we need to include this to the duplex indentification method
                                if int(prompt_pos)>=overlap_start and int(prompt_pos)<=overlap_end:
                                    #flag='YES'
                                    #here it means that we consider this pair of reads
                                    readPos=0
                                    refPos=int(POS)
                                    for eachValue,eachFlag in zip(numValues,flagValues):
                                        #print('%s with %s '%(eachValue,eachFlag)),
                                        insertionEscape=0
                                        for idx in range(1,int(eachValue)+1):
                                            if eachFlag=='S':
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t-\t%s\t-\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,eachFlag,SEQ[readPos],QUAL[readPos],extraField))
                                                readPos=readPos+1
                                            elif eachFlag=='M':
                                                position=refPos
                                                refBase=myFasta.fetch(RNAME,position-1,position)
                                                varFlag=''
                                                if refBase.upper()==SEQ[readPos].upper():
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                                else:
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                                readPos=readPos+1
                                                refPos=refPos+1
                                            elif eachFlag=='I':
                                                if insertionEscape==0:
                                                    myRefPos=int(refPos)-1
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t%s\t%s\t-\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,str(myRefPos),eachFlag,SEQ[readPos],QUAL[readPos],extraField))
                                                    readPos=readPos+1
                                                    insertionEscape=insertionEscape+1
                                                else:
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t-\t%s\t-\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,eachFlag,SEQ[readPos],QUAL[readPos],extraField))
                                                    readPos=readPos+1
                                                    insertionEscape=insertionEscape+1
                                            elif eachFlag=='D':
                                                position=refPos
                                                refBase=myFasta.fetch(RNAME,position-1,position)
                                                OutFile.write('%s\t%s\t%s\t%s\t-\t%d\t%s\t%s\t-\tNO\t-\t%s\n'%(posKey,QNAME,RNAME,MAPQ,refPos,eachFlag,refBase.upper(),extraField))
                                                refPos=refPos+1
                                            elif eachFlag=='N':
                                                position=refPos
                                                refBase=myFasta.fetch(RNAME,position-1,position)
                                                OutFile.write('%s\t%s\t%s\t%s\t-\t%d\t%s\t%s\t-\tNO\t-\t%s\n'%(posKey,QNAME,RNAME,MAPQ,refPos,eachFlag,refBase.upper(),extraField))
                                                refPos=refPos+1
                                            elif eachFlag=='=':
                                                position=refPos
                                                refBase=myFasta.fetch(RNAME,position-1,position)
                                                varFlag=''
                                                if refBase.upper()==SEQ[readPos].upper():
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                                else:
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                                readPos=readPos+1
                                            elif eachFlag=='X':
                                                position=refPos
                                                refBase=myFasta.fetch(RNAME,position-1,position)
                                                varFlag=''
                                                if refBase.upper()==SEQ[readPos].upper():
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                                else:
                                                    OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                                readPos=readPos+1
                                            elif eachFlag=='H':
                                                countH=countH+1
                                            elif eachFlag=='P':
                                                countP=countP+1
                                    OutFile.write('\n')
                                #used for debbuging
                                #if countH>0 or countP>0:
                                    #print('\t\tWarning: Found %d positions with CIGAR flag H and %d positions with CIGAR flag P.'%(countH,countP))
                                #OutFile.close()
                                #else:
                                    #flag='NO'

                                #print('QNAME:%s\tVariantPos:%s\t [%d - %d] and [%d - %d] with [%d - %d] %s\n'%(QNAME,prompt_pos,a,b,c,d,overlap_start,overlap_end,flag))
                                #print('QNAME:%s\tCIGAR:%s\tSEQ:%s\tREAD_LEN:%d\tnumValues:%s\tflagValues:%s\tINS_LEN:%s --> %d and %s'%(QNAME,CIGAR,SEQ,len(SEQ),numValues,flagValues,TLEN_int,prompt_pos,eachValue))
                            else:
                                #same concept as before
                                readPos=0
                                refPos=int(POS)
                                for eachValue,eachFlag in zip(numValues,flagValues):
                                    #print('%s with %s '%(eachValue,eachFlag)),
                                    insertionEscape=0
                                    for idx in range(1,int(eachValue)+1):
                                        if eachFlag=='S':
                                            OutFile.write('%s\t%s\t%s\t%s\t%d\t-\t%s\t-\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,eachFlag,SEQ[readPos],QUAL[readPos],extraField))
                                            readPos=readPos+1
                                        elif eachFlag=='M':
                                            position=refPos
                                            refBase=myFasta.fetch(RNAME,position-1,position)
                                            varFlag=''
                                            if refBase.upper()==SEQ[readPos].upper():
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                            else:
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                            readPos=readPos+1
                                            refPos=refPos+1
                                        elif eachFlag=='I':
                                            if insertionEscape==0:
                                                myRefPos=int(refPos)-1
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t%s\t%s\t-\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,str(myRefPos),eachFlag,SEQ[readPos],QUAL[readPos],extraField))
                                                readPos=readPos+1
                                                insertionEscape=insertionEscape+1
                                            else:
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t-\t%s\t-\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,eachFlag,SEQ[readPos],QUAL[readPos],extraField))
                                                readPos=readPos+1
                                                insertionEscape=insertionEscape+1
                                        elif eachFlag=='D':
                                            position=refPos
                                            refBase=myFasta.fetch(RNAME,position-1,position)
                                            OutFile.write('%s\t%s\t%s\t%s\t-\t%d\t%s\t%s\t-\tNO\t-\t%s\n'%(posKey,QNAME,RNAME,MAPQ,refPos,eachFlag,refBase.upper(),extraField))
                                            refPos=refPos+1
                                        elif eachFlag=='N':
                                            position=refPos
                                            refBase=myFasta.fetch(RNAME,position-1,position)
                                            OutFile.write('%s\t%s\t%s\t%s\t-\t%d\t%s\t%s\t-\tNO\t-\t%s\n'%(posKey,QNAME,RNAME,MAPQ,refPos,eachFlag,refBase.upper(),extraField))
                                            refPos=refPos+1
                                        elif eachFlag=='=':
                                            position=refPos
                                            refBase=myFasta.fetch(RNAME,position-1,position)
                                            varFlag=''
                                            if refBase.upper()==SEQ[readPos].upper():
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                            else:
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                            readPos=readPos+1
                                        elif eachFlag=='X':
                                            position=refPos
                                            refBase=myFasta.fetch(RNAME,position-1,position)
                                            varFlag=''
                                            if refBase.upper()==SEQ[readPos].upper():
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tNO\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                            else:
                                                OutFile.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\tYES\t%s\t%s\n'%(posKey,QNAME,RNAME,MAPQ,readPos,refPos,eachFlag,refBase.upper(),SEQ[readPos],QUAL[readPos],extraField))
                                            readPos=readPos+1
                                        elif eachFlag=='H':
                                            countH=countH+1
                                        elif eachFlag=='P':
                                            countP=countP+1
                                OutFile.write('\n')
                            break
            if countH>0 or countP>0:
                skipPrint=0
                #print('\t\tWarning: Found %d positions with CIGAR flag H and %d positions with CIGAR flag P.'%(countH,countP))
            OutFile.close()
        else:
            ts = time.time()
            st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('\n\t\t[%s] ERROR from function generateFlatFile: The reference genome file does not exist!\n'%(st))
            print('************************************************************************************************************************************\n')
            sys.exit()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function generateFlatFile: The input sam file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()

#parse the position list and generate one file per position with the reads that support the variants
def splitFlatFile(threadName,posKey,flatFileName,SAM_FILES):

    if exitFlag:
      threadName.exit()
    #check if the flatFileName is valid
    if os.path.exists(flatFileName):
        #the file exists so parse the positions one by one
        #print('\t\tProcessing position:%s\n'%(posKey))
        command='grep -e '+posKey+' '+flatFileName+' > '+SAM_FILES+'/'+posKey+'_reads.txt'
        os.system(command)

    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t\t[%s] ERROR from function processFlatFile: The input flat file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit() 

def variantScanner(posList,bamFilePrefix,RESULTS,SAM_FILES,index):

    #store the QNAMES of all reads with variant for the positions of interest
    #the key will be the chrom and the value will be the actual QNAME
    qnameDict=defaultdict(list)

    #prepare the output file and write it directly to out folder
    outFileName=RESULTS+'/'+bamFilePrefix+'_VariantReport.txt'
    outFile=open(outFileName,'w')
    outFile.write('#CHROM\tPOSITION\tREFERENCE\tTOTAL_READS\tA\tC\tG\tT\tDEL\tDist,Frag_A,DUP_A\tDist,Frag_C,DUP_C\tDist,Frag_G,DUP_G\tDist,Frag_T,DUP_T\tDist,Frag_DEL,DUP_DEL\tDist,Frag_INS,DUP_INS\n')
    #parse all positions
    for pkey in posList:
        tmp=pkey.split('_')
        promptChrom=tmp[0]
        promptPos=tmp[1]
        currentFile=SAM_FILES+'/'+pkey+'_reads.txt'
        if os.path.exists(currentFile):
            fileIN=open(currentFile,'r')
            readList=[]
            allReads=0
            mutatedReads=0
            promptReference=''
            variantDict=defaultdict(list)
            
            for eachLine in fileIN:
                line=eachLine.rstrip('\n')
                line=line.rstrip('\r')
                tmp=line.split('\t')
                posKey=tmp[0]
                QNAME=tmp[1]
                chrom=tmp[2]
                MAPQ=tmp[3]
                readPos=tmp[4]
                readRef=tmp[5]
                cigarFlag=tmp[6]
                refBase=tmp[7]
                currentBase=tmp[8]
                variantFlag=tmp[9]
                QUAL=tmp[10]
                duplexInfo=tmp[11]
                #print('%s\t%s\n'%(QNAME,variantFlag))
                if(chrom==promptChrom and readRef==promptPos):
                    if refBase!='-':
                        promptReference=refBase
                    if(variantFlag=='YES'):
                        allReads=allReads+1
                        mutatedReads=mutatedReads+1
                        key=currentBase
                        variantDict[key].append(eachLine)
                    elif(variantFlag=='NO' and cigarFlag=='D'):
                        #allReads=allReads+1
                        key='D'
                        variantDict[key].append(eachLine)
                    elif(variantFlag=='NO' and cigarFlag=='I'):
                        key='I'
                        variantDict[key].append(eachLine)
                    elif(variantFlag=='NO' and cigarFlag=='M'):
                        variantDict[promptReference].append(eachLine)
                        allReads=allReads+1

            readsA=variantDict['A']
            readsC=variantDict['C']
            readsG=variantDict['G']
            readsT=variantDict['T']
            readsD=variantDict['D']
            readsI=variantDict['I']
            #used to debug
            #print(readsI)
            #print('%s Total reads: %d and mutated: %d with A=%d C=%d G=%d and T=%d'%(pkey,allReads,mutatedReads,readsA,readsC,readsG,readsT))
            fileIN.close()
            #first find the reference
            if(promptReference=='A'):
                #reference is A so process C,G,T
                myStrA=processVariantReads(readsA,bamFilePrefix,RESULTS,pkey,'A')
                myStrC=processVariantReads(readsC,bamFilePrefix,RESULTS,pkey,'C')
                myStrG=processVariantReads(readsG,bamFilePrefix,RESULTS,pkey,'G')
                myStrT=processVariantReads(readsT,bamFilePrefix,RESULTS,pkey,'T')
                myStrD=processVariantReads(readsD,bamFilePrefix,RESULTS,pkey,'D')
                myStrI=processVariantReads(readsI,bamFilePrefix,RESULTS,pkey,'I')
                outFile.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n'%(promptChrom,promptPos,promptReference,allReads,len(readsA),len(readsC),len(readsG),len(readsT),len(readsD),len(readsI),myStrA,myStrC,myStrG,myStrT,myStrD,myStrI))
                for j in readsA:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsC:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsG:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsT:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsD:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsI:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)

            elif(promptReference=='C'):
                #reference is C so process A,G,T
                myStrA=processVariantReads(readsA,bamFilePrefix,RESULTS,pkey,'A')
                myStrC=processVariantReads(readsC,bamFilePrefix,RESULTS,pkey,'C')
                myStrG=processVariantReads(readsG,bamFilePrefix,RESULTS,pkey,'G')
                myStrT=processVariantReads(readsT,bamFilePrefix,RESULTS,pkey,'T')
                myStrD=processVariantReads(readsD,bamFilePrefix,RESULTS,pkey,'D')
                myStrI=processVariantReads(readsI,bamFilePrefix,RESULTS,pkey,'I')
                outFile.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n'%(promptChrom,promptPos,promptReference,allReads,len(readsA),len(readsC),len(readsG),len(readsT),len(readsD),len(readsI),myStrA,myStrC,myStrG,myStrT,myStrD,myStrI))
                for j in readsA:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsG:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsT:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsD:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsI:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)

            elif(promptReference=='G'):
                #reference is G so process A,C,T
                myStrG=processVariantReads(readsG,bamFilePrefix,RESULTS,pkey,'G')
                myStrC=processVariantReads(readsC,bamFilePrefix,RESULTS,pkey,'C')
                myStrA=processVariantReads(readsA,bamFilePrefix,RESULTS,pkey,'A')
                myStrT=processVariantReads(readsT,bamFilePrefix,RESULTS,pkey,'T')
                myStrD=processVariantReads(readsD,bamFilePrefix,RESULTS,pkey,'D')
                myStrI=processVariantReads(readsI,bamFilePrefix,RESULTS,pkey,'I')
                outFile.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n'%(promptChrom,promptPos,promptReference,allReads,len(readsA),len(readsC),len(readsG),len(readsT),len(readsD),len(readsI),myStrA,myStrC,myStrG,myStrT,myStrD,myStrI))
                for j in readsC:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsA:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsT:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsD:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsI:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)

            elif(promptReference=='T'):
                #reference is T so process A,C,T
                myStrT=processVariantReads(readsT,bamFilePrefix,RESULTS,pkey,'T')
                myStrC=processVariantReads(readsC,bamFilePrefix,RESULTS,pkey,'C')
                myStrG=processVariantReads(readsG,bamFilePrefix,RESULTS,pkey,'G')
                myStrA=processVariantReads(readsA,bamFilePrefix,RESULTS,pkey,'A')
                myStrD=processVariantReads(readsD,bamFilePrefix,RESULTS,pkey,'D')
                myStrI=processVariantReads(readsI,bamFilePrefix,RESULTS,pkey,'I')
                outFile.write('%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n'%(promptChrom,promptPos,promptReference,allReads,len(readsA),len(readsC),len(readsG),len(readsT),len(readsD),len(readsI),myStrA,myStrC,myStrG,myStrT,myStrD,myStrI))
                for j in readsC:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsG:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsA:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsD:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)
                for j in readsI:
                    line=j.rstrip('\n')
                    line=line.rstrip('\r')
                    tmp=line.split('\t')
                    Key=tmp[2]
                    Value=tmp[1]
                    qnameDict[Key].append(Value)

            del variantDict
            
        else:
            ts=time.time()
            st=datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
            print('\n\t\t[%s] ERROR from function processFlatFile: The input flat file does not exist!\n'%(st))
    outFile.close()

    #erase the _read.txt file to save some space
    command='rm '+currentFile
    #dont run
    #os.system(command)

    return qnameDict

#this functions processes the reads carrying the variant of interest and reports the information we need
def processVariantReads(myHash,bamFilePrefix,RESULTS,pkey,base):

    #variable dummy
    dummy=0
    XI=''
    DUP=0
    positionHashNoXI=defaultdict(list)
    positionHashXI=defaultdict(list)
    positionHash=defaultdict(list)

    extraPositionHash=defaultdict(list)

    for recordKey in myHash:
        line=recordKey.rstrip('\n')
        line=line.rstrip('\r')
        tmp=line.split('\t')
        QNAME=tmp[1]
        myInfo=tmp[11]
        myInfo=myInfo.split(':')
        POS=myInfo[1]
        RNEXT=myInfo[2]
        PNEXT=myInfo[3]
        TLEN=myInfo[4]
        XI=int(myInfo[5])
        FLAG=myInfo[6]
        #consider only those that have XI>2 to count them
        if(XI>=2):
            posKey=POS+'_'+str(abs(int(TLEN)))
            positionHashXI[posKey].append(QNAME)

            R1=int(POS)
            R2=int(PNEXT)

            firstValue=min(R1,R2)
            secondValue=max(R1,R2)
            posKey=str(firstValue)+'_'+str(secondValue)+'_'+str(abs(int(TLEN)))
            extraPositionHash[posKey].append(QNAME)

        else:
            posKey=POS+'_'+str(abs(int(TLEN)))
            positionHashNoXI[posKey].append(QNAME)
            #print('%s %s '%(pkey,QNAME))   
        #store the starting positions of all reads, the TLEN and the FLAG to process the duplexes
        if RNEXT=='=':
            R1=int(POS)
            R2=int(PNEXT)
            if R1<R2:
                positionKey=POS
                positionValue=TLEN+'_'+FLAG
                positionHash[positionKey].append(positionValue)
            elif R1==R2:
                dummy=0
            else:
                #here we have situation where we see the mate so we store as leftmost position the PNEXT
                #as TLEN the abs(TLEN) and the inverse FLAG
                positionKey=PNEXT
                myLen=abs(int(TLEN))
                currentFlag=int(FLAG)
                b=64
                c=128
                myFlag=''
                if ((currentFlag&b)==64 and (currentFlag&c==0)):
                #the current flag is first in pair to the other read is the second
                    myFlag=128
                elif (currentFlag&b)==0 and (currentFlag&c==128):
                #the current flag is second in pair so make it first
                    myFlag=64
                positionValue=str(myLen)+'_'+str(myFlag)
                positionHash[positionKey].append(positionValue)
        #at this point we have all info we need to find duplexes
    for dictIdx in positionHash:
        recordsFound=positionHash[dictIdx]
        count=len(recordsFound)
        #compute the duplexes per position, need to index them by TLEN
        insHash=defaultdict(list)
        for TlenIdx in recordsFound:
            tmp=TlenIdx.split("_")
            TLEN=tmp[0]
            FLAG=tmp[1]
            insKey=dictIdx+'_'+TLEN
            insHash[insKey].append(TlenIdx)
        #dupl=findPairOrientation(insHash,dictOUT)
        dupl=findPairOrientation(insHash)
        DUP=DUP+dupl
        del insHash
    #dictOUT.close()
    #print(extraPositionHash)
    myStr=str(len(positionHashXI))+','+str(len(extraPositionHash))+','+str(DUP)
    return myStr

def findPairOrientation(insHash):
    #scan the Hash and find all duplexes per start position
    DUP=0
    for insSize in insHash:
        recordsFound=insHash[insSize]
        firstCount=0
        secondCount=0
        #duplexes per TLEN
        duplexes=0
        for idx1 in recordsFound:
            #scan the reads and find how many are first in pair and how many are second in pair
            #print('%s '%idx1),
            tmp=idx1.split("_")
            TLEN=tmp[0]
            FLAG=tmp[1]
            a=int(FLAG)
            b=64
            c=128
            #print(' *%s'%FLAG),
            if ((a&b)==64 and (a&c==0)):
                firstCount=firstCount+1
            elif (a&b)==0 and (a&c==128):
                secondCount=secondCount+1
        duplexes=min(firstCount,secondCount)
        DUP=DUP+duplexes
        #dictOUT.write('%s\t%s\t%s\n'%(insSize,len(recordsFound),duplexes))
    return DUP

#erase the big sam files to save some space
def cleanFiles(SAM_FILES):
    command='rm '+SAM_FILES+'/*.txt && rm '+SAM_FILES+'/*.bed && rm '+SAM_FILES+'/*.sam && rm '+SAM_FILES+'/*.bam'
    os.system(command)

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=6:
        printUsage()
    else:
        bamFile=sys.argv[1].split('bamFile=')
        bamFile=bamFile[1]
        
        bedFile=sys.argv[2].split('positionFile=')
        bedFile=bedFile[1]
        
        refGenome=sys.argv[3].split('referenceGenome=')
        refGenome=refGenome[1]
        #read the output dir
        outDIR=sys.argv[4].split('outDIR=')
        outDIR=outDIR[1]
        
        index=sys.argv[5].split('index=')
        index=index[1]

        bamFilePrefix=bamFile[:-4]
        tmp=bamFilePrefix.split('/')
        bamFilePrefix=tmp[-1]
        
        #save the bam file prefix for further naming of files
        #generate folders for storing the final results and the intermediate results
        RESULTS=outDIR+'/'+bamFilePrefix+'/RESULTS_duplexCaller'
        SAM_FILES=outDIR+'/'+bamFilePrefix+'/INTERM_FILES_duplexCaller'
        bamFilePrefix_DIR=outDIR+'/'+bamFilePrefix
        command='mkdir -p '+bamFilePrefix_DIR+' && mkdir -p '+RESULTS+' && mkdir -p '+SAM_FILES
        os.system(command)

        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        #print('\n\t[%s] Processing sample %s'%(st,bamFilePrefix))

        #store the positions for testing
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Function storePositionsFile: store the input genomic positions'%(st))
        posList=storePositionsFile(bedFile)
        chromList=storeChromosomes(bedFile)

        #splitting of chromosomes
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Spliting the original BAM per chromosome: Releasing threads...'%(st))
        threads = []
        i=1
        for chrom in chromList:
        # Create new threads and start them
            myStr='Thread-'+str(i)
            thread = myThread(i, myStr,chrom,bamFilePrefix,bamFile,SAM_FILES)
            thread.start()
            threads.append(thread)
            i=i+1
        #find your threads
        for t in threads:
            t.join()
        del threads

        #collect them
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Spliting the original BAM per chromosome: All threads collected...'%(st))

        #generate a sam file with all reads spanning the position of interest
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        targetSamFileName=SAM_FILES+'/'+bamFilePrefix+'_target_f'+index+'.mysam'
        print('\n\t[%s] Function generateTargetSamFile...'%(st))
        generateTargetSamFile(posList,bamFile,targetSamFileName,SAM_FILES,index)

        flatFileName=SAM_FILES+'/'+bamFilePrefix+'_TargetFlat_f'+index+'.txt'
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Function generateFlatFile...'%(st))
        generateFlatFile(targetSamFileName,flatFileName,refGenome)

        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Produce list of reads supporting the input variants: Releasing threads...'%(st))
        #release the threads
        threadLock = threading.Lock()
        threads = []
        i=1
        for posKey in posList:
        # Create new threads and start them
            myStr='Thread-'+str(i)
            thread = myVariantThread(i, myStr,posKey,flatFileName,SAM_FILES)
            thread.start()
            threads.append(thread)
            i=i+1
        #find threads
        for t in threads:
            t.join()
        #collect them
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Produce list of reads supporting the input variants: All threads collected...'%(st))

        #parse the read files and generate variant report
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Function variantScanner: identify the variants'%(st))
        QNAME_DICT=variantScanner(posList,bamFilePrefix,outDIR,SAM_FILES,index)

        #clear the intermediate files
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n\t[%s] Clear all intermediate files produced.'%(st))
        cleanFiles(SAM_FILES)


#start 
if __name__=='__main__':
    myMain()


