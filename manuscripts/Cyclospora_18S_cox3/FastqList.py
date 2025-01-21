#!/usr/bin/env python

import argparse,os
import pandas as pd
from functools import partial, reduce

##### Create a fastq list of representative sequences that meet a community threshold (percent as decimal)
##### Community thresholds determined based on duplicate count details (generated from seqkit rmdup)
##### qyr9@cdc.gov

def getFiles(directory, suffix='_duplicated.detail.txt'):
    """Create list of files and sample names from seqkit outputs"""
    fileList = []
    nameList = []
    for f in os.listdir(os.path.abspath(directory)):
        if (f.endswith(suffix)):
            fileList.append(f)
            fns = f.split(suffix)[0]
            fns2 = str(fns).replace("-","_")
            nameList.append(fns2)
    return fileList, nameList

def makedict(fileList,nameList,directory):
    """Create dictionary of duplicate results"""
    dictcounts = {}
    for i in range(len(fileList)):
        dictcounts[nameList[i]] = pd.read_csv(os.path.abspath(directory) +'/'+ fileList[i], sep='\t', header=None,
                names=['counts','seqs'])
    return dictcounts 

def getthresholds(dictcounts, percentage):
    """Create a dictionary of reads above a % community threshold"""
    dictreads = {}
    dictreadcounts = {}
    
    for key, df in dictcounts.items():
        df['seqs'] = df['seqs'].str.replace(r"\,.*","", regex=True)
        total = df['counts'].sum()
        threshold = round(total*percentage)
        newdf = df[df['counts'] >= threshold]
        dictreadcounts[key] = newdf
        seqlist = newdf[['seqs']]
        dictreads[key] = seqlist
        
    return dictreads, dictreadcounts

def getoutput(dictreads,dictreadcounts=""):
    """Output sequence df from dictionary"""
    
    for key, df in dictreads.items(): 
        df.to_csv(key+'_fastq-list.txt', sep='\t', index=False, header=False)
   
    if not dictreadcounts:
        pass
    else:
        for key, df in dictreadcounts.items(): 
            df.to_csv(key+'_counts.txt', sep='\t', index=False, header=False) 

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Create a fastq list of seqs above community threshold', epilog='_____')
    parser.add_argument('-d','--directory',type=str,required=True,help="directory containing files")
    parser.add_argument('-s','--suffix',type=str,required=True,help="file name suffix for duplicated details")
    parser.add_argument('-p','--percentage',type=float,required=True,default='0',help= "threshold percentage as decimal")

    args = parser.parse_args()


    if args.directory == '.':
        directory = os.getcwd()
    else:
        directory = args.directory


    fileList, nameList = getFiles(directory, args.suffix)
    dictcounts = makedict(fileList, nameList,directory)
    dictreads, dictreadcounts = getthresholds(dictcounts,args.percentage)
    getoutput(dictreads,dictreadcounts)
