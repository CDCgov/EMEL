#!/usr/bin/env python

import os, argparse
import pandas as pd
from functools import partial, reduce

##### 07/2023
##### Parse blast output and returns total counts (unique top hits + randomized tied top hits based on max bitscores)
##### Requires blast ouptut outfmt 6 with the following categories: qseqid sacc staxid stitle length pident qcovs mismatch gapopen qstart qend sstart send evalue bitscore
##### Outputs based on customized blast databases: Isolate (Salmonella LT2), ATCC Environmental, and Zymo
##### qyr9@cdc.gov

def getFiles(directory, suffix1):
    """Create list of files and sample names from blast results"""
    fileList = []
    nameList = []
    for f in os.listdir(os.path.abspath(directory)):
        if (f.endswith(suffix1)):
            fileList.append(f)
            fns = f.split(suffix1)[0]
            nameList.append(fns)
    return fileList, nameList

def makedict(fileList,nameList,directory):
    """Convert blast results to dictionary of dfs, requires blast results with following categories listed below"""
    blastdict = {}
    for i in range(len(fileList)):
        blastdict[nameList[i]] = pd.read_csv(os.path.abspath(directory) + '/' + fileList[i], sep='\t', header=None,
                names=['seqID','acc','taxID','title','length','PID','qcov','mismatch','gapopen','qstart',
                      'qend','sstart','send','evalue','bitscore'])
    return blastdict       

def addtaxa(blastdict,community=""):
    """ Create identical taxID for taxa hits from same genome from blast search"""
    taxdict = {}

    if community == 'Environmental':
        for key, df in blastdict.items():
            df.loc[df['acc'].str.contains('Staph'), 'taxID'] = 'S.epidermis'
            df.loc[df['acc'].str.contains('Salmonella'), 'taxID'] = 'Salmonella'
            df.loc[df['acc'].str.contains('coli'), 'taxID'] = 'Ecoli'
            df.loc[df['acc'].str.contains('Chromo'), 'taxID'] = 'C.violaceum'
            df.loc[df['acc'].str.contains('haloplankt'), 'taxID'] = 'P.haloplanktis'
            df.loc[df['acc'].str.contains('subtilis'), 'taxID'] = 'B.subtilis'
            df.loc[df['acc'].str.contains('fluor'), 'taxID'] = 'P.fluorescens'
            df.loc[df['acc'].str.contains('faecalis'), 'taxID'] = 'E.faecalis'
            df.loc[df['acc'].str.contains('halophilus'), 'taxID'] = 'H.halophilus'
            df.loc[df['acc'].str.contains('luteus'), 'taxID'] = 'M.luteus'

            df.loc[df['acc'].str.contains('uv'), 'taxID'] = 'Vector'
            df.loc[df['acc'].str.contains('Bacillus_phage'), 'taxID'] = 'Bacillus_phage'
            df.loc[df['acc'].str.contains('phi29'), 'taxID'] = 'Phi29'
            df.loc[df['acc'].str.contains('phiX'), 'taxID'] = 'PhiX'
            taxdict[key] = df
    
        return taxdict
    
    if community == 'Zymo':
        for key, df in blastdict.items():
            df.loc[df['acc'].str.contains('subtilis'), 'taxID'] = 'B.subtilis'
            df.loc[df['acc'].str.contains('Crypto'), 'taxID'] = 'C.neoformans'
            df.loc[df['acc'].str.contains('faecalis'), 'taxID'] = 'E.faecalis'
            df.loc[df['acc'].str.contains('coli'), 'taxID'] = 'Ecoli'
            df.loc[df['acc'].str.contains('fermentum'), 'taxID'] = 'L.fermentum'
            df.loc[df['acc'].str.contains('Listeria'), 'taxID'] = 'L.monocytogenes'
            df.loc[df['acc'].str.contains('aeruginosa'), 'taxID'] = 'P.aeruginosa'
            df.loc[df['acc'].str.contains('Saccharo'), 'taxID'] = 'S.cerevisiae'
            df.loc[df['acc'].str.contains('Salmonella'), 'taxID'] = 'Salmonella'
            df.loc[df['acc'].str.contains('aureus'), 'taxID'] = 'S.aureus'

            df.loc[df['acc'].str.contains('uv'), 'taxID'] = 'Vector'
            df.loc[df['acc'].str.contains('Bacillus_phage'), 'taxID'] = 'Bacillus_phage'
            df.loc[df['acc'].str.contains('phi29'), 'taxID'] = 'Phi29'
            df.loc[df['acc'].str.contains('phiX'), 'taxID'] = 'PhiX'
            
        return taxdict
    
    if community == 'Isolate':
        for key, df in blastdict.items():
            df.loc[df['acc'].str.contains('coli'), 'taxID'] = 'Ecoli'
            df.loc[df['acc'].str.contains('Salmonella'), 'taxID'] = 'Salmonella'
            df.loc[df['acc'].str.contains('sapiens'), 'taxID'] = 'Human'

            df.loc[df['acc'].str.contains('uv'), 'taxID'] = 'Vector'
            df.loc[df['acc'].str.contains('Bacillus_phage'), 'taxID'] = 'Bacillus_phage'
            df.loc[df['acc'].str.contains('phi29'), 'taxID'] = 'Phi29'
            df.loc[df['acc'].str.contains('phiX'), 'taxID'] = 'PhiX'
            
        return taxdict
    
    else:   
        print('Not a listed community')

def gethits(taxdict):
    """Create dictionary of unique and tied hits based on max bit score"""
    maxscoredict = {}
    uniquehitsdict = {}
    tiedhitsdict = {}
    tiedtaxadict = {}
    
    val = 1
 
    for key, df in blastdict.items():
        scores = df[df.bitscore == df.bitscore.groupby(df['seqID']).transform('max')]
        maxscoredict[key] = scores
        
    for key, df in maxscoredict.items():      
        nodups = df.drop_duplicates(['seqID','taxID','bitscore'], keep='first')
        unique = nodups.drop_duplicates(subset=['seqID','bitscore'], keep=False)
        uniquehits = unique.assign(count=[val for w in range(len(unique.index))])
        uniquehitsdict[key] = uniquehits
        
        ties = df.duplicated(['seqID','bitscore'], keep=False)
        tiedf = df[ties].sort_values('seqID').reset_index(drop=True)
        ties2 = tiedf.drop_duplicates(['seqID','taxID','bitscore'], keep='first')
        resort = ties2.sort_values(['seqID', 'taxID'])
        tiedhitsdict[key] = resort
    
    for key, df in tiedhitsdict.items():
        taxa_series = df.groupby('seqID')['taxID'].apply('+'.join)
        taxa_df = taxa_series.to_frame().reset_index()
        taxa_counts = taxa_df.assign(count=[val for w in range(len(taxa_df.index))])
        tiedtaxadict[key] = taxa_counts
        
    return uniquehitsdict, tiedhitsdict, tiedtaxadict
    
def getcounts(hitsdict):
    """ Return counts dictionary"""
    finalcountsdict = {}
    
    for key, df in hitsdict.items():
        taxid = df[['taxID','count']]
        totals = taxid.groupby(['taxID'],as_index=False).sum(numeric_only=True)
        totalsrename = totals.rename(columns = {'title':'Taxa','count':key})
        finalcountsdict[key] = totalsrename
        
    return(finalcountsdict)

def taxonomy(countsdict):
    """Convert dictionaries to df and return full taxonomy and reduced taxonomy dfs"""
    
    df_reduce = partial(pd.merge, on = ['taxID'], how='outer')
    tax = reduce(df_reduce,countsdict.values())
    tax.fillna(0, inplace=True)

    return(tax)

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Parse BLAST directory and return tophits based on max bitscore', epilog='_____')
    parser.add_argument('-d','--directory',type=str,required=True,help="directory containing BLAST output files")
    parser.add_argument('-s','--suffix',type=str,required=True,help="file name suffix, use -s='-suffix' if suffix starts with '-'")
    parser.add_argument('-mock', '--community',required=True,default="",help="mock community database used in search")

    args = parser.parse_args()

    if args.directory == '.':
        directory = os.getcwd()
    else:
        directory = args.directory

    fileList, fileName = getFiles(directory, args.suffix)
    blastdict = makedict(fileList, fileName, directory)
    taxdict = addtaxa(blastdict, args.community)
    uniquehitsdict, tiedhitsdict, tiedtaxadict = gethits(taxdict)

    uniquecountsdict = getcounts(uniquehitsdict)
    tiedcountsdict = getcounts(tiedtaxadict)

    uniquetax = taxonomy(uniquecountsdict)
    tiedtax = taxonomy(tiedcountsdict)
        
    uniquetax.to_csv("Unique_taxa_hits.tsv", sep='\t', index=False)
    tiedtax.to_csv("Tied_taxa_hits.tsv", sep='\t', index=False)