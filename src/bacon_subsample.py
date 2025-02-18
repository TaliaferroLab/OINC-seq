#As a statistical framework for identifying genes with differing 8OG-mediated conversion rates
#across conditions, use a subsampling approach.  For each gene, subsample the reads assigned to it,
#calculating a porc value for each subsample.  So then you end up with a distribution of porc
#values for each gene in each sample.

#This relies on two pickled dictionaries produced by pigpen.py: read2gene.pkl and readconvs.pkl.
#The first is a dictionary of the form {readid : ensembl_gene_id}.
#The second is a dictionary of the form {readid : {convs}}
#where {convs} is of the form {x_y : count} where x is the reference sequence and y is the query sequence.
#x is one of 'agct' and y is one of 'agctn'.

#Then, using Hoteling t-test, compare distributions of porc values across conditions.

import pickle
import numpy as np
import random
import sys
from collections import defaultdict

def getReadsperGene(read2genepkl):
    #pigpen produces a dictionary of the form {readid : ensembl_gene_id} (assignreads.processOverlaps)
    #It has written this dictionary to a pickled file.
    #We need a dictionary of the form {ensembl_gene_id : [readids]}

    print('Loading gene/read assignments...')
    with open(read2genepkl, 'rb') as infh:
        read2gene = pickle.load(infh)
    print('Done!')

    readspergene = {} #{ensembl_gene_id : [readids that belong to this gene]}

    readcount = 0
    for read in read2gene:
        gene = read2gene[read]
        if gene not in readspergene:
            readspergene[gene] = [read]
        else:
            readspergene[gene].append(read)

    return readspergene

def makeGeneConvdict(readspergene, readconvspkl):
    #Want to make a dictionary that looks like this
    #{gene : [{convs1}, {convs2}, ...]} where each conv dict corresponds to one read
    #that has been assigned to this gene

    print('Loading read conversion info...')
    with open(readconvspkl, 'rb') as infh:
        readconvs = pickle.load(infh)
    print('Done!')

    print('Making gene : read conversion dictionary...')
    geneconv = defaultdict(list)
    for gene in readspergene:
        reads = readspergene[gene]
        for read in reads:
            try:
                geneconv[gene].append(readconvs[read])
            except KeyError:
                pass #a read for which we calculated conversions but didn't get assigned to any read
    print('Done!')

    return geneconv

def calcPORC(convs):
    #accepting a list of conversion dictionaries

    convGcount = 0
    totalGcount = 0
    allconvcount = 0
    allnonconvcount = 0

    for conv in convs:
        convG = conv['g_t'] + conv['g_c']
        totalG = conv['g_t'] + conv['g_c'] + \
            conv['g_a'] + conv['g_n'] + conv['g_g']

        allconv = conv['a_t'] + conv['a_c'] + conv['a_g'] + conv['g_t'] + conv['g_c'] + \
            conv['g_a'] + conv['t_a'] + conv['t_c'] + \
            conv['t_g'] + conv['c_t'] + conv['c_g'] + conv['c_a'] + \
            conv['a_n'] + conv['g_n'] + conv['c_n'] + conv['t_n']

        allnonconv = conv['a_a'] + conv['g_g'] + conv['c_c'] + conv['t_t']

        convGcount += convG
        totalGcount += totalG
        allconvcount += allconv
        allnonconvcount += allnonconv

    allnt = allconvcount + allnonconvcount
    try:
        convGrate = convGcount / totalGcount
    except ZeroDivisionError:
        convGrate = np.nan

    try:
        totalmutrate = allconvcount / allnt
    except ZeroDivisionError:
        totalmutrate = np.nan

    #Calculate porc
    if totalmutrate == np.nan:
        porc = np.nan
    elif totalmutrate > 0:
        try:
            porc = np.log2(convGrate / totalmutrate)
        except:
            porc = np.nan
    else:
        porc = np.nan

    return porc

def subsamplegeneconv(geneconv, subsamplesize, n_subsamples):
    subsampledporcs = defaultdict(list)  # {ensembl_gene_id : [porc values]}

    for gene in geneconv:
        print(gene)
        convs = geneconv[gene]
        nreads = len(convs)
        n_readstosubsample = int(nreads * subsamplesize)
        for i in range(n_subsamples):
            subsampledconvs = random.sample(convs, n_readstosubsample)
            porc = calcPORC(subsampledconvs)
            print(porc)
            subsampledporcs[gene].append(porc)

    return subsampledporcs

#Take in a sampconds, calculate subsamples for each, end up with a dictionary like so:
#{gene : {condition : [[subsampled porcs sample 1], [subsampled porcs sample 2], ...]}}




if __name__ == '__main__':

    readspergene = getReadsperGene(sys.argv[1])
    geneconv = makeGeneConvdict(readspergene, sys.argv[2])
    subsamplegeneconv(geneconv, 0.3, 100)

