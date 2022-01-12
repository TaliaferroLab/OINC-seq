#Once we have determined the number of conversions in each read (getmismatches.py)
#and the gene to which each read belongs (assignreads.py),
#combine the two to get per-gene conversion rates
import sys
import pickle
from collections import Counter
import numpy as np


def getPerGene(convs, reads2gene):
    #convs = {} #{readid : {a_a : count, a_t : count, etc.}}
    #reads2gene = {} #{readid : ensembl_gene_id}

    #Reorient reads2gene so that we can get the number of reads assigned to each gene
    genes = list(reads2gene.values())
    numreadspergene = dict(Counter(genes)) #{ensembl_gene_id : numreads}

    convsPerGene = {} # {ensembl_gene_id : {conversion : count}}
    
    #Initialize convsPerGene
    geneids = reads2gene.values()
    geneids = list(set(geneids))
    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n']

    #It's possible (but relatively rare) for a read to be in convs but
    #not in reads2gene (or vice versa). Filter for reads only present in both.
    convreads = set(convs.keys())
    genereads = set(reads2gene.keys())
    commonreads = convreads.intersection(genereads)
    print('Calculating conversion rates in {0} reads.'.format(len(commonreads)))

    convs = {key:value for (key, value) in convs.items() if key in commonreads}
    reads2gene = {key:value for (key, value) in reads2gene.items() if key in commonreads}


    for gene in geneids:
        convsPerGene[gene] = {}
        for conv in possibleconvs:
            convsPerGene[gene][conv] = 0

    for read in reads2gene:
        gene = reads2gene[read]
        readConvs = convs[read]
        for c in readConvs:
            convCount = readConvs[c]
            convsPerGene[gene][c] += convCount

    return numreadspergene, convsPerGene

def writeConvsPerGene(numreadspergene, convsPerGene, outfile, use_g_t, use_g_c):
    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n']

    with open(outfile, 'w') as outfh:
        #total G is number of ref Gs encountered
        #convG is g_t + g_c (the ones we are interested in)
        outfh.write(('\t').join(['Gene', 'numreads'] + possibleconvs + ['totalG', 'convG', 'convGrate', 'G_Trate', 'G_Crate', 'porc']) + '\n')
        genes = sorted(convsPerGene.keys())

        for gene in genes:
            numreads = numreadspergene[gene]
            convcounts = []
            c = convsPerGene[gene]
            for conv in possibleconvs:
                convcount = c[conv]
                convcounts.append(convcount)

            convcounts = [str(x) for x in convcounts]

            totalG = c['g_g'] + c['g_c'] + c['g_t'] + c['g_a'] + c['g_n']
            if use_g_t and use_g_c:
                convG = c['g_c'] + c['g_t']
            elif use_g_c and not use_g_t:
                convG = c['g_c']
            elif use_g_t and not use_g_c:
                convG = c['g_t']
            elif not use_g_t and not use_g_c:
                print('ERROR: we have to be counting either G->T or G->C, if not both!')
                sys.exit()
                
            g_ccount = c['g_c']
            g_tcount = c['g_t']

            totalmut = c['a_t'] + c['a_c'] + c['a_g'] + c['g_t'] + c['g_c'] + c['g_a'] + c['t_a'] + c['t_c'] + c['t_g'] + c['c_t'] + c['c_g'] + c['c_a']
            totalnonmut = c['a_a'] + c['g_g'] + c['c_c'] + c['t_t']
            allnt = totalmut + totalnonmut

            try:
                convGrate = convG / totalG
            except ZeroDivisionError:
                convGrate = 'NA'
                
            try: 
                g_crate = g_ccount / totalG
            except ZeroDivisionError:
                g_crate = 'NA'

            try:
                g_trate = g_tcount / totalG
            except ZeroDivisionError:
                g_trate = 'NA'

            try:
                totalmutrate = totalmut / allnt
            except ZeroDivisionError:
                totalmutrate = 'NA'

            #normalize convGrate to rate of all mutations
            #Proportion Of Relevant Conversions
            if totalmutrate == 'NA':
                porc = 'NA'
            elif totalmutrate > 0:
                try:
                    porc = np.log2(convGrate / totalmutrate)
                except:
                    porc = 'NA'
            else:
                porc = 'NA'

            #Format numbers for printing
            if type(convGrate) == float:
                convGrate = '{:.2e}'.format(convGrate)
            if type(g_trate) == float:
                g_trate = '{:.2e}'.format(g_trate)
            if type(g_crate) == float:
                g_crate = '{:.2e}'.format(g_crate)
            if type(porc) == np.float64:
                porc = '{:.3f}'.format(porc)

            outfh.write(('\t').join([gene, str(numreads)] + convcounts + [str(totalG), str(convG), str(convGrate), str(g_trate), str(g_crate), str(porc)]) + '\n')

        
if __name__ == '__main__':
    with open('OINC3.mDBF.subsampled.filtered.convs.pkl', 'rb') as infh:
        convs = pickle.load(infh)
    with open('read2gene.pkl', 'rb') as infh:
        reads2gene = pickle.load(infh)
    numreadspergene, convsPerGene = getPerGene(convs, reads2gene)
    writeConvsPerGene(numreadspergene, convsPerGene, 'test.txt')


    