#!/usr/bin/env python

import os,argparse
import pandas as pd
import numpy as np
from functools import partial, reduce

##### Parse blast output and returns total counts (unique top hits + randomized tied top hits based on max bitscores)
##### Requires blast ouptut outfmt 6 with the following categories: qseqid sacc staxid stitle length pident qcovs mismatch gapopen qstart qend sstart send evalue bitscore
##### Requires unique fasta reads blast outputs and duplicated.detail outputs from seqkit rmdup; each sample must have a duplicated.detail output (if no duplicates in sample, create dummy duplicate detail txt file)
##### Optional taxa target search (returns all unique taxaID hits + all taxaID hits tied with other taxa). 
##### qyr9@cdc.gov

def getFiles(directory, suffix1, suffix2):
    """Create list of files and sample names from blast results and rep dupcounts from gene slices"""
    fileList1 = []
    nameList1 = []
    fileList2 = []
    nameList2 = []
    
    for f in os.listdir(os.path.abspath(directory)):
        if (f.endswith(suffix1)):
            fileList1.append(f)
            fns = f.split(suffix1)[0]
            nameList1.append(fns)
    
    for f in os.listdir(os.path.abspath(directory)):
        if (f.endswith(suffix2)):
            fileList2.append(f)
            fns = f.split(suffix2)[0]
            nameList2.append(fns)
            
    return fileList1, nameList1, fileList2, nameList2

def makedict(fileList1,nameList1,fileList2,nameList2,directory):
    """Convert data files to dataframe dictionaries"""
    blastdict = {}
    dupsdict = {}
    
    for i in range(len(fileList1)):
        blastdict[nameList1[i]] = pd.read_csv(os.path.abspath(directory) +'/'+ fileList1[i], sep='\t', header=None,
                names=['seqID','acc','taxID','title','length','PID','qcov','mismatch','gapopen','qstart',
                      'qend','sstart','send','evalue','bitscore'])
        
    for i in range(len(fileList2)):
        dupsdict[nameList2[i]] = pd.read_csv(os.path.abspath(directory) +'/'+ fileList2[i], sep='\t', header=None,
                names=['count','seqID'])
        
    return blastdict, dupsdict 

def gethits(blastdict, target=""):
    """Create dictionary of unique, tied and combined top hits based on max bit score"""
    maxscoredict = {}
    uniquehitsdict = {}
    tiedhitsdict = {}
    tophitsdict = {}
    tiedhitsdictall = {}
    
    val = 1
 
    for key, df in blastdict.items():
        scores = df[df.bitscore == df.bitscore.groupby(df['seqID']).transform('max')]
        maxscoredict[key] = scores
        
    if not target:
        for key, df in maxscoredict.items():
            nodups = df.drop_duplicates(['seqID','bitscore'], keep=False)
            unique = nodups.drop_duplicates(subset=['seqID'], keep='first')
            uniquehitsdict[key] = unique
        
            ties = df.duplicated(['seqID','bitscore'], keep=False)
            dfrand = df[ties].sample(frac=1).reset_index(drop=True)
            tiedhitsall = dfrand.sort_values('seqID').reset_index(drop=True)
            tiedhitsdictall[key] = tiedhitsall
        
            tied = tiedhitsall.drop_duplicates(subset=['seqID'], keep='first')
            tiedhitsdict[key] = tied
        
            totalhits = pd.concat([unique,tied], axis=0)
            tophitsdict[key] = totalhits
        
        return uniquehitsdict, tiedhitsdictall, tiedhitsdict, tophitsdict
    
    else:
        
        if target.isdigit():
            targetID = int(target)
            
            for key, df in maxscoredict.items():
                
                targetrows = df[df['taxID'] == targetID]
                totalhits = targetrows.drop_duplicates(['seqID','bitscore'], keep='first')
                tophitsdict[key] = totalhits
            
                notaxadups = df.drop_duplicates(['seqID','taxID','bitscore'], keep='first') 
                target_unique = notaxadups.drop_duplicates(subset=['seqID','bitscore'], keep=False)
                unique = target_unique[target_unique['taxID'] == targetID]
                uniquehitsdict[key] = unique
            
                dups = df.duplicated(['seqID','bitscore'], keep=False)
                dups_sort = df[dups].sort_values('seqID').reset_index(drop=True)
                taxa_tied = dups_sort.drop_duplicates(['seqID','taxID','bitscore'], keep='first')
                one_taxa = taxa_tied[taxa_tied['taxID'].isin([targetID]).groupby(taxa_tied['seqID']).transform('any')]
            
                dups2 = one_taxa.duplicated(['seqID','bitscore'], keep= False)
                tiedhitsall = one_taxa[dups2].sort_values('seqID').reset_index(drop=True)
                tiedhitsdictall[key] = tiedhitsall
                
                duprows = tiedhitsall[tiedhitsall['taxID'] == targetID]
                tied = duprows.drop_duplicates(['seqID','bitscore'], keep=False)
                tiedhitsdict[key] = tied
                
        else:
            
            for key, df in maxscoredict.items():
            
                targetrows = df[df.title.str.contains(target)]
                totalhits = targetrows.drop_duplicates(['seqID','taxID','bitscore'], keep='first')
                tophitsdict[key] = totalhits
            
                notaxadups = df.drop_duplicates(['seqID','taxID','bitscore'], keep='first') 
                target_unique = notaxadups.drop_duplicates(subset=['seqID','bitscore'], keep=False)
                unique = target_unique[target_unique.title.str.contains(target)]
                uniquehitsdict[key] = unique
            
                dups = df.duplicated(['seqID','bitscore'], keep= False)
                dups_sort = df[dups].sort_values('seqID').reset_index(drop=True)
                taxa_tied = dups_sort.drop_duplicates(['seqID','taxID','bitscore'], keep='first') 
                one_taxa = taxa_tied[taxa_tied['title'].str.contains(target).groupby(taxa_tied['seqID']).transform('any')]
            
                dups2 = one_taxa.duplicated(['seqID','bitscore'], keep= False)
                tiedhitsall = one_taxa[dups2].sort_values('seqID').reset_index(drop=True)
                tiedhitsdictall[key] = tiedhitsall
                
                duprows = tiedhitsall[tiedhitsall.title.str.contains(target)]
                tiedhitsdict[key] = duprows
        
        return uniquehitsdict, tiedhitsdictall, tiedhitsdict, tophitsdict

def modifydetails(dupsdict):
    """ Remove duplicate sequence IDs from dupcount detail"""
    dupsdictmod = {}
    
    for key, df in dupsdict.items():
        df['seqID'] = df['seqID'].str.replace(r"\,.*","", regex=True)
        dupsdictmod[key] = df
        
    return dupsdictmod
    
def mergedict(hitsdict,dupsdictmod,combine='taxID'):
    """Merge BLAST hits dictionary with duplicate counts dictionary"""
    countsdict = {}
    tiedcountsdict = {}
    
    hitsdict = {str(i):dict_i for i, dict_i in enumerate([hitsdict,dupsdictmod])} 
    hitsdictcomb = {group_key:[df_dict[group_key] for df_dict in hitsdict.values()] 
                      for group_key in list(hitsdict.values())[0].keys()}
    hitsdictmerge = {group_key:reduce(lambda left,right:
                                pd.merge(left,right, on=['seqID'], how='left'),
                                dat_df_list)
                       for group_key,dat_df_list in hitsdictcomb.items()}
     
    for key, df in hitsdictmerge.items():
        replace = df['count'].fillna(1, inplace=True)
    
    if 'taxID' in combine:
        for key, df in hitsdictmerge.items(): 
            if not df.empty:
                taxid = df[['taxID','title','count']]
                totals = taxid.groupby(['taxID','title'],as_index=False).sum(numeric_only=False)
                totalsrename = totals.rename(columns = {'title':'Taxa','count':key})
                countsdict[key] = totalsrename
    else:
        for key, df in hitsdictmerge.items():
            if not df.empty:
                seqid = df[['seqID','taxID','title','count']]
                countsdict[key] = seqid

    return countsdict

def taxonomy(countsdict,toptax=""):
    """Convert dictionaries to df and return full taxonomy and reduced taxonomy dfs"""
    if not countsdict:
        dummydf = pd.DataFrame(0, index=np.arange(1), columns=np.arange(len(toptax.columns)))
        dummydf.columns = toptax.columns
        return(dummydf,dummydf)
    else:
        df_reduce = partial(pd.merge, on = ['taxID','Taxa'], how='outer')
        tax = reduce(df_reduce,countsdict.values())
        tax.fillna(0, inplace=True)
        taxsum = tax.groupby(['taxID']).sum(numeric_only=True)
        taxid = tax[['taxID','Taxa']]
        nodups = taxid.drop_duplicates(subset=['taxID'], keep='first')
        reducetax =  pd.merge(nodups,taxsum, on=['taxID'], how='outer')
        
        return(tax,reducetax)

def sortcolumns(tax,reducetax):
    """Sort columns by R1 and R2 if present"""
    if tax.columns.str.contains('R1').any():
        fullinfo = tax[['taxID','Taxa']]
        fullsort = tax.loc[:, ~tax.columns.isin(["taxID","Taxa"])].sort_index(key=lambda x: x.str.split('R1').str[1], axis=1)
        fullfinal = pd.concat([fullinfo,fullsort], axis=1)
        reduceinfo = reducetax[['taxID','Taxa']]
        reducesort = reducetax.loc[:, ~reducetax.columns.isin(["taxID","Taxa"])].sort_index(key=lambda x: x.str.split('R1').str[1], axis=1)
        reducefinal = pd.concat([reduceinfo,reducesort], axis = 1)
        return(fullfinal,reducefinal)
    else:
        return(tax,reducetax)

def combinedf(topsortreduce,uniquesortreduce,tiedsortreduce):
    """Combined dataframes"""
    topsortreduce['label'] ='all hits (total)'
    uniquesortreduce['label'] = 'unique hits (by taxID)'
    tiedsortreduce['label'] = 'tied hits (with other taxID)'
    Taxadf = topsortreduce[['taxID','Taxa']]
    
    df1 = topsortreduce.drop(['Taxa'], axis = 1)
    df2 = uniquesortreduce.drop(['Taxa'], axis=1)
    df3 = tiedsortreduce.drop(['Taxa'], axis = 1)
    
    newdf = pd.concat([df1,df2,df3])
    newdf.fillna(0, inplace=True)
    mergedf = pd.merge(newdf, Taxadf, how="left", on=["taxID"])
    order = ['taxID','Taxa','label']
    finaldf = mergedf[order + [c for c in mergedf if c not in order]]
    sortdf = finaldf.sort_values(by=['taxID','label'])

    return(sortdf)

def dicttodf(countsdict):
    """Output raw tiedcount df from dictionary"""
    for key, df in countsdict.items(): 
        df.to_csv(key +'_top_tied_hits_all.txt', sep='\t', index=False, header=False)


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Parse BLAST directory and return tophits based on max bitscore', epilog='_____')
    parser.add_argument('-d','--directory',type=str,required=True,help="directory containing BLAST output files")
    parser.add_argument('-s1','--suffix1',type=str,required=True,help="file name suffix for blast files, use -s='-suffix' if suffix starts with '-'")
    parser.add_argument('-s2','--suffix2',type=str,required=True,help="file name suffix for duplicate count gene slice files, use -s='-suffix' if suffix starts with '-'")
    parser.add_argument('-tie', '--tiedhits',required=False,action='store_true',help="returns all top tied hit dfs for all samples and combined count table")
    parser.add_argument('-unique', '--uniquehits',required=False,action='store_true',help="returns unique tophits combined count table")
    parser.add_argument('-sort', '--sortcolumns',required=False,action='store_true',help="sorts count tables by R1 and R2")
    parser.add_argument('-target', '--target',required=False,default="",help="taxaID or taxa string to search")
    parser.add_argument('-o', '--output_prefix', required=False,default="",help="output prefix")

    args = parser.parse_args()

    if args.directory == '.':
        directory = os.getcwd()
    else:
        directory = args.directory

    fileList1, fileName1, fileList2, fileName2 = getFiles(directory, args.suffix1, args.suffix2)
    blastdict, dupsdict  = makedict(fileList1, fileName1, fileList2, fileName2, directory)
    uniquehitsdict, tiedhitsdictall, tiedhitsdict, tophitsdict = gethits(blastdict,args.target)
    dupsdictmod = modifydetails(dupsdict)
    
    topcountsdict = mergedict(tophitsdict,dupsdictmod,'taxID')
    toptax, topreducetax = taxonomy(topcountsdict)
    topsort, topsortreduce = sortcolumns(toptax, topreducetax)


    if args.target == "":
        topsort.to_csv(args.output_prefix + "tophits_full_tax.tsv", sep='\t', index=False)
        topsortreduce.to_csv(args.output_prefix + "tophits_reduce_tax.tsv", sep='\t', index=False)
    
        if args.tiedhits:
            tiedcountsdict = mergedict(tiedhitsdict,dupsdictmod,'taxID')
            tiedcountsdictall = mergedict(tiedhitsdictall,dupsdictmod,'seqID')
            tiedtax, tiedreducetax = taxonomy(tiedcountsdict,toptax)
            tiedsort, tiedsortreduce = sortcolumns(tiedtax, tiedreducetax)
            tiedsort.to_csv(args.output_prefix +"tiedhits_full_tax.tsv", sep='\t', index=False)
            tiedsortreduce.to_csv(args.output_prefix + "tiedhits_reduce_tax.tsv", sep='\t', index=False)
            dicttodf(tiedhitsdictall)
        
        else:
            pass
    
        if args.uniquehits:
            uniquecountsdict = mergedict(uniquehitsdict,dupsdictmod,'taxID')
            uniquetax, uniquereducetax = taxonomy(uniquecountsdict,toptax)
            uniquesort, uniquesortreduce = sortcolumns(uniquetax, uniquereducetax)
            uniquesort.to_csv(args.output_prefix + "uniquehits_full_tax.tsv", sep='\t', index=False)
            uniquesortreduce.to_csv(args.output_prefix + "uniquehits_reduce_tax.tsv", sep='\t', index=False)

        else:
            pass

    else:
        #topsortreduce.to_csv(args.output_prefix + "alltophits.tsv", sep='\t', index=False)

        uniquecountsdict = mergedict(uniquehitsdict,dupsdictmod,'taxID')
        uniquetax, uniquereducetax = taxonomy(uniquecountsdict,toptax)
        uniquesort, uniquesortreduce = sortcolumns(uniquetax, uniquereducetax)

        tiedcountsdict = mergedict(tiedhitsdict,dupsdictmod,'taxID')
        tiedcountsdictall = mergedict(tiedhitsdictall,dupsdictmod,'seqID')
        tiedtax, tiedreducetax = taxonomy(tiedcountsdict,toptax)
        tiedsort, tiedsortreduce = sortcolumns(tiedtax, tiedreducetax)

        combined = combinedf(topsortreduce,uniquesortreduce,tiedsortreduce)
        combined.to_csv(args.output_prefix + "tophits_combined.tsv", sep='\t', index = False)

        if args.tiedhits:
            for key, df in tiedhitsdictall.items(): 
                df.to_csv(key + '_' + args.output_prefix + 'tied_taxa.tsv', sep='\t', index=False, header=False)
            

        else:
            pass

        if args.uniquehits:
            uniquetax.to_csv(args.output_prefix + "unique_taxa.tsv", sep='\t', index=False)

        else:
            pass
            
        
        
