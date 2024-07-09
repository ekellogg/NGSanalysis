
from Bio.Seq import Seq

import dataclasses
import pandas as pd
import gzip
import matplotlib as plt
from collections import Counter
import numpy as np
from joblib import Parallel, delayed
import time
import gc
from functools import reduce
import operator
import itertools
import multiprocessing
from multiprocessing import Manager

def MutationInfo(seq1, seq2, delimiter="+"):
    differences = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
    out = ""
    if len(differences) == 0:
        return "WT"
    for dd in differences:
        out = out + seq1[int(dd)] + str(dd+1) + seq2[int(dd)]
        if len(differences) > 1:
            out = out + delimiter

    if len(differences) > 1:
        return out[:-1]
    return out


def findOffset(primer, seq):
    # Find the starting index of the last occurrence of substring1 in string2
    start_index = seq.rfind(primer)
    
    if start_index == -1:
        # substring1 is not found in string2
        return -1
    
    # Compute the index of the last character of substring1 within string2
    last_char_index = start_index + len(primer)
    return last_char_index


def NumberMutations(seq1,seq2):
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def translateDNA(ntSeq):	
    dna_seq = Seq(ntSeq)
    return str(dna_seq.translate())

def flattenList(lst):
    return reduce(operator.add, lst)

#templateSeq needs to be nucleotide, not protein
def readNGSfiles(file, templateSeq, fwdPrimer,nworkers=12):
    
    with open(file) as f:
        lines = f.readlines()
    manager = Manager()
    reads = manager.list(lines[1::4]) 
    # Find the offset of the primer in the first sequence line
    offset = findOffset(fwdPrimer, lines[1])
    if offset == -1:
        raise ValueError("Forward primer not found in the sequence")
    refAAseq = translateDNA(templateSeq)

    print('number of lines: ' + str(len(lines)) + ' and reads: ' + str(len(reads)) + '\n')
    del lines
    gc.collect()
    
    if nworkers==0:
        results = readLines(reads,templateSeq,offset,0,len(reads),file)    
    else:
        nitems = int(len(reads)/(nworkers))
        StartNdx = list(range(0,len(reads),nitems))
        startTime = time.time()
        results = Parallel(nworkers)(delayed(readLines)(reads,templateSeq,offset,indices,nitems,file) for indices in StartNdx)
        endTime = time.time()
        print('time to process file: '+str(endTime-startTime))

    variants = []
    protseq = []
    filename = []
    fileline = []
    NTmutinfo = []
    numNTmuts = []
    AAmutinfo = []
    numAAmuts = []
    
    if nworkers == 0:
        variants = results[0]
        protseq = results[1]
        filename = results[2]
        fileline = results[3]
        NTmutinfo = results[4]
        numNTmuts = results[5]
        AAmutinfo = results[6]
        numAAmuts = results[7]
    else:

        for i in range(0,nworkers):
             variants = [*variants, *results[i][0]]
             protseq = [*protseq, *results[i][1]]
             filename = [*filename, *results[i][2]]
             fileline = [*fileline, *results[i][3]]
             NTmutinfo = [*NTmutinfo, *results[i][4]]
             numNTmuts = [*numNTmuts, *results[i][5]]
             AAmutinfo = [*AAmutinfo, *results[i][6]]
             numAAmuts = [*numAAmuts, *results[i][7]]
    
    df = {
          'variants': pd.Series(variants),
          'protseq': pd.Series(protseq),
          'filename': pd.Series(filename),       
          'fileline': pd.Series(fileline),
          'NTmutinfo': pd.Series(NTmutinfo),
          'numNTmuts': pd.Series(numNTmuts),
          'AAmutinfo': pd.Series(AAmutinfo),
          'numAAmuts': pd.Series(numAAmuts)
    }

    return pd.DataFrame(df)
	

#if using parallel, should use existing protein sequences
def readLines(lines,templateSeq,offset,indexStart,nitems,filelabel):

    readlen = len(templateSeq)  
    variants = []
    protseq = []
    fileName = []
    fileLine = []
    NTmutinfo = []
    numNTmuts = []
    AAmutinfo = []
    numAAmuts = []
    refAAseq = translateDNA(templateSeq)
    endndx = indexStart+nitems
    if endndx > len(lines):
        endndx = len(lines)
    for i in range(indexStart, endndx):  # Get sequences and counts
        seq_i = lines[i].rstrip()[offset:(offset + readlen)]
        if len(seq_i) == 498: #check that read is not truncated

            aaSeq = translateDNA(seq_i)
            AA = MutationInfo(aaSeq,refAAseq)
            NT = MutationInfo(seq_i,templateSeq)
            variants.append(seq_i)
            protseq.append(aaSeq)
            fileName.append(filelabel)
            fileLine.append(i)
            NTmutinfo.append(NT)
            AAmutinfo.append(AA)
            numNTmuts.append(NumberMutations(seq_i,templateSeq))
            numAAmuts.append(NumberMutations(aaSeq,refAAseq))

    return [variants,protseq,fileName,fileLine,NTmutinfo,numNTmuts,AAmutinfo,numAAmuts]
  

def filterByCounts(df,minvalue = 100,nucleotide=True):
    c = None
    if nucleotide:
        c = df['variants'].value_counts()
        df['counts'] = df['variants'].map(c) #adds a new column called 'counts'
    else:
        c = df['protseq'].value_counts()
        df['counts'] = df['protseq'].map(c)
        
    return df[df['counts'] > minvalue]

def removeWT(df,templateSeq,nucleotide=True):
    if nucleotide:
        labels = df.NTmutinfo
        return df[labels != "WT"]
    else:
        labels = df.AAmutinfo
        return df[labels != "WT"]

def getCounts(df,nucleotide=True):
    if nucleotide:
        c = df['variants'].value_counts()
        df['counts'] = df['variants'].map(c) #adds a new column called 'counts'
    else:
        c = df['protseq'].value_counts()
        df['counts'] = df['protseq'].map(c)
    return df

def OneHotEncoding(sequence, amino_acids='ACDEFGHIKLMNPQRSTVWYX*'):

    seqlen = len(sequence)
    
    # Create a matrix of zeros with shape (len(sequence), len(amino_acids))
    one_hot_matrix = np.zeros((seqlen, len(amino_acids)), dtype='int')
    
    # Create a dictionary to map each amino acid to an index
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}
    
    # Fill in the one-hot matrix        
    for i, aa in enumerate(sequence):
        if aa in aa_to_index:
            one_hot_matrix[i, aa_to_index[aa]] =+ 1
        else:
            raise ValueError(f"Unknown amino acid: {aa}, from {sequence}")
    
    return one_hot_matrix

def findCommonVariants(df1,df2,nucleotide=True):
    k = set()
    if nucleotide:
        k = set(df1['variants']).intersection(set(df2['variants']))
        return [k,df1[df1['variants'].isin(k)], df2[df2['variants'].isin(k)]]
    else:
        k = set(df1['protseq']).intersection(set(df2['protseq']))
        return [df1[df1['protseq'].isin(k)], df2[df2['protseq'].isin(k)]]
    
    
# Function to compute coverage for a column
def compute_coverage(df1,df2,col):
    unique_values_df1 = set(df1[col].unique())
    unique_values_df2 = set(df2[col].unique())
    common_values = unique_values_df1.intersection(unique_values_df2)
    coverage = len(common_values) / len(unique_values_df1)
    return coverage

#compute log2 enrichments from dataframes
def getEnrichments(df1,df2,mincounts = 100,nucleotide=True):
    col = 'protseq'
    if nucleotide:
        col = 'variants'
        
    df1filt = filterByCounts(df1,mincounts,nucleotide)
    df2filt = filterByCounts(df2,mincounts,nucleotide)
    df1counts = df1filt[col].value_counts()
    df2counts = df2filt[col].value_counts()
    aligned1, aligned2 = df1counts.align(df2counts, 'outer',fill_value=1e-6)
    relative_frequencies = aligned2 / aligned1
    finite_mask = np.isfinite(relative_frequencies)
    filtered_freq = relative_frequencies[finite_mask]
    enrichments= np.log2(filtered_freq)
    return enrichments
    
def removeNan(df):
    mask = np.isfinite(df) & df.notna() 
    filtered_df = df[mask.all(axis=1)]
    filtered_df = filtered_df[ (filtered_df > -20).all(axis=1) ] 
    filtered_df = filtered_df[ (filtered_df < 20).all(axis=1) ]
    return filtered_df

