#!/usr/local/bin/python

'''
				Parse a variant report returned by duplexCaller and filter the results to produce a readable report
                      
BEGIN COPYRIGHT NOTICE

    cfDNApipeline code -- (c) 2019 Dimitrios Kleftogiannis -- GIS -- A*STAR
   	
   	Copyright 2019 Genome Institute of Singapore (GIS) and Agency for Science, Technology and Research (A*STAR).

    This Program is free software licensed under the MIT License.

    You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENSE.md".

    This Program is distributed as a service to the research community and is experimental in nature and may have hazardous properties. 
    
    The Program is distributed WITHOUT ANY WARRANTY, express or implied. In particular all warranties as to SATISFACTORY QUALITY or FITNESS FOR A PARTICULAR PURPOSE are excluded. 

    Published reports of research using this code (or a modified version) should cite the relevant article of this tool.

    Comments and bug reports are welcome.
       
    Email to dimitrios.kleftogiannis@gkaust.edu.sa

    I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
    You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE

'''

import sys
import os
import re
from collections import defaultdict
from itertools import groupby
import datetime
import time

#prints information about program's execution
def printUsage():
    print('To run this program please type the following:')
    print('\tpython filterVariantReportAdjusted.py inputFile=report.txt originalVCF=vcf.txt\n')
    print('Where:\n') 
    print('\treport.txt is a variant report returned by duplexCaller\n')
    print('\tvcf.txt is the original VCF file that contains reference and mutated allele\n')
    print('Please give the arguments in the indicated order similar to the provided example!\n') 


def storeVariantReport(inputFile,originalVCF,inputFilePrefix):

	idx=0
	variantDict=defaultdict(list)
	if os.path.exists(inputFile):
	#the file exists
		InFile=open(inputFile,'r')
		for eachLine in InFile:
			if idx==0:
				line = eachLine.rstrip('\n')
				tmp=line.split("\t")
				idx=idx+1
			else:
				line = eachLine.rstrip('\n')
				tmp=line.split("\t")
				chrom = tmp[0]
				pos = tmp[1]
				refBase = tmp[2]
				totalReads = tmp[3]
				readsA = tmp[4]
				readsC = tmp[5]
				readsG = tmp[6]
				readsT = tmp[7]
				readsDEL = tmp[8]
				readsINS=tmp[9]
				infoA = tmp[10]
				infoC = tmp[11]
				infoG = tmp[12]
				infoT = tmp[13]
				infoDEL = tmp[14]
				infoINS = tmp[15]

				#store per position each nucleotide
				
				infoA = infoA.split(',')
				distFragA = infoA[0]
				distPairsA = infoA[1]
				duplexA = infoA[2]
				keyA=chrom+'_'+pos+'_A'
				AF_A=round(float(readsA)/float(totalReads),5)
				valueA=totalReads+'\t'+readsA+'\t'+str(AF_A)+'\t'+distFragA+'\t'+distPairsA+'\t'+duplexA
				variantDict[keyA].append(valueA)

				infoC = infoC.split(',')
				distFragC = infoC[0]
				distPairsC = infoC[1]
				duplexC = infoC[2]
				keyC=chrom+'_'+pos+'_C'
				AF_C=round(float(readsC)/float(totalReads),5)
				valueC=totalReads+'\t'+readsC+'\t'+str(AF_C)+'\t'+distFragC+'\t'+distPairsC+'\t'+duplexC
				variantDict[keyC].append(valueC)

				infoG = infoG.split(',')
				distFragG = infoG[0]
				distPairsG = infoG[1]
				duplexG = infoG[2]
				keyG=chrom+'_'+pos+'_G'
				AF_G=round(float(readsG)/float(totalReads),5)
				valueG=totalReads+'\t'+readsG+'\t'+str(AF_G)+'\t'+distFragG+'\t'+distPairsG+'\t'+duplexG
				variantDict[keyG].append(valueG)

				infoT = infoT.split(',')
				distFragT = infoT[0]
				distPairsT = infoT[1]
				duplexT = infoT[2]
				keyT=chrom+'_'+pos+'_T'
				AF_T=round(float(readsT)/float(totalReads),5)
				valueT=totalReads+'\t'+readsT+'\t'+str(AF_T)+'\t'+distFragT+'\t'+distPairsT+'\t'+duplexT
				variantDict[keyT].append(valueT)

				#info for deletions
				infoDEL = infoDEL.split(',')
				distFragDEL = infoDEL[0]
				distPairsDEL = infoDEL[1]
				duplexDEL = infoDEL[2]
				keyDEL=chrom+'_'+pos+'_DEL'
				AF_DEL=round(float(readsDEL)/(float(readsDEL)+float(totalReads)),5)
				#AF_DEL=round(float(readsDEL)/float(totalReads),5)
				valueDEL=totalReads+'\t'+readsDEL+'\t'+str(AF_DEL)+'\t'+distFragDEL+'\t'+distPairsDEL+'\t'+duplexDEL
				variantDict[keyDEL].append(valueDEL)

				#info for insertions
				infoINS = infoINS.split(',')
				distFragINS = infoINS[0]
				distPairsINS = infoINS[1]
				duplexINS = infoINS[2]
				keyINS=chrom+'_'+pos+'_INS'
				AF_INS=round(float(readsINS)/(float(readsINS)+float(totalReads)),5)
				#AF_INS=round(float(readsINS)/float(totalReads),5)
				valueINS=totalReads+'\t'+readsINS+'\t'+str(AF_INS)+'\t'+distFragINS+'\t'+distPairsINS+'\t'+duplexINS
				variantDict[keyINS].append(valueINS)

		InFile.close()

		#we have the positions of interest stored and scan the original file to report the results
	if os.path.exists(originalVCF):
		InFile=open(originalVCF,'r')
		idx=0
		outFileName=inputFilePrefix+'_VariantsClean.txt'
		outFile=open(outFileName,'w')
		for eachLine in InFile:
			if idx==0:
				line = eachLine.rstrip('\n')
				line=line.rstrip()
				idx=idx+1
				outFile.write('%s\n'%(line))
			elif idx==1:
				#line = eachLine.rstrip('\n')
				#tmp=line.split("\t")
				idx=idx+1
				outFile.write('#CHROM\tPOS\tREF\tALT\tGENE\tvariantType\tCOV\tSupporintReads\tVAF\tRefBias\tVarBias\tMeanReadPos\tMeanReadQual\tMeanMappingQual\tSignalToNoise\tHighQualVarReads\tHighQualCov\tConsequence\tImpact\tBiotype\tProteinPos\tAminoAcidChange\tExistingVar\tPopulationVAF\tPredictedClinicalSignif\tMyTotalReads\tMySupportingReads\tMyAF\tDistReads\tDistReadPairs\tDuplexes\tAdjustedVAF\tRefDuplexes\tPASS\n')
			else:
				idx=idx+1
				line = eachLine.rstrip('\n')
				tmp=line.split("\t")
				
				#read the first parts
				CHROM = tmp[0]
				POS = tmp[1]
				REF = tmp[2]
				ALT = tmp[3]
				GENE = tmp[4]
				variantType = tmp[5]
				dp = int(tmp[6])
				vd = int(tmp[7])
				vaf = float(tmp[8])

				refbias = tmp[9]
				varbias = tmp[10]
				t=varbias.split(':')
				FW=int(t[0])
				BW=int(t[1])
				pmean = float(tmp[11])
				qual = float(tmp[12])
				mq = float(tmp[13])
				sn = float(tmp[14])
				hicnt = int(tmp[15])
				hicov = int(tmp[16])

				Consequence = tmp[17]
				IMPACT = tmp[18]
				BIOTYPE = tmp[19]
				Protein_position = tmp[20]
				Amino_acids = tmp[21]
				Existing_variation = tmp[22]
				PopulationVAF = tmp[23]
				CLIN_SIG = tmp[24]

				flag='tio'

				if variantType=='SNV':
					myPos=pos
					promptKey=CHROM+'_'+POS+'_'+ALT
					recordFound=variantDict[promptKey]
					#check if found
					if len(recordFound)!=0:
						for myIDX in recordFound:
							tmp=myIDX.split('\t')
							totalReads=int(tmp[0])
							supportReads=int(tmp[1])
							AF=float(tmp[2])
							distFrag=int(tmp[3])
							distPairs=int(tmp[4])
							duplexes=int(tmp[5])

							A=0
							V=0
							K=0
							N=0
							a=0

							referencePromptkey=CHROM+'_'+POS+'_'+REF
							refRecordFound = variantDict[referencePromptkey]
							for my in refRecordFound:
								tmp=my.split('\t')
								refTotalReads=tmp[0]
								refSupportReads=tmp[1]
								refAF=tmp[2]
								refDistFrag=int(tmp[3])
								refdistPairs=int(tmp[4])
								refDuplexes=int(tmp[5])

								A=duplexes
								V=int(supportReads)
								K=refDuplexes
								N=int(refSupportReads)

							if (V+N-A-K)!=0:
								AdjustedVAF=float(V-A)/(float(V+N-A-K))
							else:
								a=0

							if (hicov>=100 and hicnt>=3 and FW>1 and BW>1 and pmean>15 and sn>20 and duplexes>=1 and distFrag>=2 and distPairs>1):
								flag='YES'
							else:
								flag='NO'
							#MyTotalReads\tMySupportingReads\tMyAF\tDistReads\tDistReadPairs\tDuplexes\tAdjustedVAF\tRefDuplexes\tPASS\n')
							outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%.3f\t%d\t%s\n'%(CHROM,POS,REF,ALT,GENE,variantType,dp,vd,vaf,refbias,varbias,pmean,qual,mq,sn,hicnt,hicov,Consequence,IMPACT,BIOTYPE,Protein_position,Amino_acids,Existing_variation,PopulationVAF,CLIN_SIG,totalReads,supportReads,AF,distFrag,distPairs,duplexes,AdjustedVAF,refDuplexes,flag))
							

					else:
						print('\n\tNot found%s\n'%(promptKey))

				elif variantType=='Deletion':
					myExtraCount=0
					c=len(REF)
					if c>20:
						print('\n\tBig deletions are skipped....\n')
						continue
					i=0
					for myCount in range(int(POS),int(POS)+c):
						myPos=myCount
						promptKey=CHROM+'_'+str(myPos)+'_DEL'
						recordFound=variantDict[promptKey]

						referencePromptkey=CHROM+'_'+str(myPos)+'_'+REF[i]
						#print('%s - %s\n'%(referencePromptkey,refBase[i]))
						
						refRecordFound = variantDict[referencePromptkey]

						if len(recordFound)!=0:
							#print('%s with %s : '%(pos,promptKey)),
							#print(recordFound)
							for myIDX in recordFound:
								tmp=myIDX.split('\t')
								totalReads=int(tmp[0])
								supportReads=int(tmp[1])
								AF=float(tmp[2])
								distFrag=int(tmp[3])
								distPairs=int(tmp[4])
								duplexes=int(tmp[5])

								A=0
								V=0
								K=0
								N=0
								a=0

								for my in refRecordFound:
									tmp=my.split('\t')
									refTotalReads=tmp[0]
									refSupportReads=tmp[1]
									refAF=tmp[2]
									refDistFrag=int(tmp[3])
									refdistPairs=int(tmp[4])
									refDuplexes=int(tmp[5])
									#a=refSupportReads+'_'+str(refDuplexes)

								A=duplexes
								V=int(supportReads)
								K=refDuplexes
								N=int(refSupportReads)
								
								if (V+N-A-K)!=0:
									AdjustedVAF=float(V-A)/(float(V+N-A-K))
								else:
									AdjustedVAF=0
								if (hicov>=200 and hicnt>=6 and FW>1 and BW>1 and pmean>15 and sn>25 and duplexes>2 and distFrag>=2 and distPairs>2):
									flag='YES'
									myExtraCount=myExtraCount+1
								else:
									flag='NO'
								ALT='-'
								ref=REF[i]
								outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%.3f\t%d\t%s\n'%(CHROM,POS,ref,ALT,GENE,variantType,dp,vd,vaf,refbias,varbias,pmean,qual,mq,sn,hicnt,hicov,Consequence,IMPACT,BIOTYPE,Protein_position,Amino_acids,Existing_variation,PopulationVAF,CLIN_SIG,totalReads,supportReads,AF,distFrag,distPairs,duplexes,AdjustedVAF,refDuplexes,flag))
						i=i+1
						if myExtraCount>0:
							flag='YES'
						else:
							flag='NO'
				elif variantType=='Insertion':
					myPos=POS
					promptKey=CHROM+'_'+myPos+'_INS'
					recordFound=variantDict[promptKey]

					referencePromptkey=CHROM+'_'+str(myPos)+'_'+REF
					refRecordFound = variantDict[referencePromptkey]
					if len(recordFound)!=0:
							#print('%s with %s : '%(pos,promptKey)),
							#print(recordFound)
							for myIDX in recordFound:
								tmp=myIDX.split('\t')
								totalReads=int(tmp[0])
								supportReads=int(tmp[1])
								AF=float(tmp[2])
								distFrag=int(tmp[3])
								distPairs=int(tmp[4])
								duplexes=int(tmp[5])

								A=0
								V=0
								K=0
								N=0
								a=0

								for my in refRecordFound:
									tmp=my.split('\t')
									refTotalReads=tmp[0]
									refSupportReads=tmp[1]
									refAF=tmp[2]
									refDistFrag=int(tmp[3])
									refdistPairs=int(tmp[4])
									refDuplexes=int(tmp[5])
									#a=refSupportReads+'_'+str(refDuplexes)

								A=duplexes
								V=int(supportReads)
								K=refDuplexes
								N=int(refSupportReads)
								
								if (V+N-A-K)!=0:
									AdjustedVAF=float(V-A)/(float(V+N-A-K))
								else:
									AdjustedVAF=0
								if (hicov>=200 and hicnt>=6 and FW>1 and BW>1 and pmean>15 and sn>25 and duplexes>2 and distFrag>=2 and distPairs>2):
									flag='YES'
								else:
									flag='NO'
			
								outFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%.3f\t%d\t%s\n'%(CHROM,POS,REF,ALT,GENE,variantType,dp,vd,vaf,refbias,varbias,pmean,qual,mq,sn,hicnt,hicov,Consequence,IMPACT,BIOTYPE,Protein_position,Amino_acids,Existing_variation,PopulationVAF,CLIN_SIG,totalReads,supportReads,AF,distFrag,distPairs,duplexes,AdjustedVAF,refDuplexes,flag))


		outFile.close()
		InFile.close()

#main function of the program
def myMain():
#check the number of input arguments
	if len(sys.argv)!=3:
		printUsage()
	else:
		 		
		inputFile=sys.argv[1].split('inputFile=')
		inputFile=inputFile[1]
		inputFilePrefix=inputFile[:-18]
		originalVCF=sys.argv[2].split('originalVCF=')
		originalVCF=originalVCF[1]

		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
		print('\n\t[%s] Function storeVariantReport for sample: %s'%(st,inputFile))
		storeVariantReport(inputFile,originalVCF,inputFilePrefix)


#start ...
if __name__=='__main__':
	myMain()




