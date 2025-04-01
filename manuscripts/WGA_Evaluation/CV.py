#!/usr/bin/env python

import argparse,os,sys
import pandas as pd
import numpy as np
from functools import partial, reduce

#### 07/2023 AA: Returns coefficient of variation of genome, chromosome, or contig 
#### based on read depth output files from samtools alignments
#### qyr9@cdc.gov

def getFiles(directory, suffix1):
    """Create list of files and sample names"""
    fileList = []
    nameList = []
    for f in os.listdir(os.path.abspath(directory)):
        if (f.endswith(suffix1)):
            fileList.append(f)
            fns = f.split(suffix1)[0]
            nameList.append(fns)
    return fileList, nameList

def makedict(fileList,nameList,directory):
    """Convert samtools depth files to dictionary"""
    dictdepth = {}
    for i in range(len(fileList)):
        dictdepth[nameList[i]] = pd.read_csv(os.path.abspath(directory) + '/' + fileList[i], sep='\t', header=None,
                names=['REF','POS','Depth'])
    return dictdepth

def getgenome(dictdepth,reference = ""):
    """Subset depth file by genome, gene, or plasmid and calculate CV of read depth"""
    genomedict = {}
    CV = {}
    cv = lambda x: np.std(x, ddof=1) / np.mean(x)
    
    for key, df in dictdepth.items():
        if not reference:
            genomedf = df
            genomedict[key] = genomedf
        else:
            genomedf = df[df.REF == reference]
            genomedict[key] = genomedf
        
    for key, df in genomedict.items():
        CoVar = cv(df["Depth"])
        CV[key] = CoVar

    CVdf = pd.DataFrame([CV],columns=CV.keys())
    outputdf = CVdf.T.reset_index()
    outputdf.columns.values[0:3] = ["Sample", "CV"]
    return outputdf


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Get coefficient of variation for genome read depth', epilog='_____')
    parser.add_argument('-d','--directory',type=str,required=True,help="directory containing samtools depth files")
    parser.add_argument('-s','--suffix',type=str,required=True,help="file name suffix")
    parser.add_argument('-r','--reference',type=str,required=False,default="",help= "Genome, gene, plasmid ID from depth file")
    #parser.add_argument('-o','--output',type=str,required=True,help= "Desired Output")

    args = parser.parse_args()

    if args.directory == '.':
        directory = os.getcwd()
    else:
        directory = args.directory
    
    fileList, fileName = getFiles(directory, args.suffix)
    dictdepth = makedict(fileList, fileName, directory)
    outputdf = getgenome(dictdepth, args.reference)

    if args.reference == "":
        outputdf.to_csv("CV.tsv", sep='\t', index=False)
    else:
        outputdf.to_csv(args.reference + "_CV.tsv", sep='\t', index=False)

