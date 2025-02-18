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
        't_a', 't_t', 't_c', 't_g', 't_n',
        'a_x', 'g_x', 'c_x', 't_x', 'ng_xg']

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

def writeConvsPerGene(sampleparams, numreadspergene, convsPerGene, outfile, use_g_t, use_g_c, use_g_x, use_ng_xg):
    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n',
        'a_x', 'g_x', 'c_x', 't_x', 'ng_xg']

    with open(outfile, 'w') as outfh:
        #Write arguments for this pigpen run
        for arg in sampleparams:
            outfh.write('#' + arg + '\t' + str(sampleparams[arg]) + '\n')
        #total G is number of ref Gs encountered
        #convG is g_t + g_c + g_x + ng_xg (the ones we are interested in)
        outfh.write(('\t').join(['Gene', 'numreads'] + possibleconvs + ['totalG', 'convG', 'convGrate', 'G_Trate', 'G_Crate', 'G_Xrate', 'NG_XGrate', 'porc']) + '\n')
        genes = sorted(convsPerGene.keys())

        for gene in genes:
            numreads = numreadspergene[gene]
            convcounts = []
            c = convsPerGene[gene]
            for conv in possibleconvs:
                convcount = c[conv]
                convcounts.append(convcount)

            convcounts = [str(x) for x in convcounts]

            totalG = c['g_g'] + c['g_c'] + c['g_t'] + c['g_a'] + c['g_n'] + c['g_x']
            convG = 0
            possiblegconv = ['g_t', 'g_c', 'g_x', 'ng_xg']
            for ind, x in enumerate([use_g_t, use_g_c, use_g_x, use_ng_xg]):
                if x == True:
                    convG += c[possiblegconv[ind]]
                
            g_ccount = c['g_c']
            g_tcount = c['g_t']
            g_xcount = c['g_x']
            ng_xgcount = c['ng_xg']

            totalmut = c['a_t'] + c['a_c'] + c['a_g'] + c['g_t'] + c['g_c'] + c['g_a'] + c['t_a'] + c['t_c'] + c['t_g'] + c['c_t'] + c['c_g'] + c['c_a'] + c['g_x'] + c['ng_xg']
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
                g_xrate = g_xcount / totalG
            except ZeroDivisionError:
                g_xrate = 'NA'

            try:
                ng_xgrate = ng_xgcount / totalG
            except ZeroDivisionError:
                ng_xgrate = 'NA'

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
            if type(g_xrate) == float:
                g_xrate = '{:.2e}'.format(g_xrate)
            if type(ng_xgrate) == float:
                ng_xgrate = '{:.2e}'.format(ng_xgrate)
            if type(porc) == np.float64:
                porc = '{:.3f}'.format(porc)

            outfh.write(('\t').join([gene, str(numreads)] + convcounts + [str(totalG), str(convG), str(convGrate), str(g_trate), str(g_crate), str(g_xrate), str(ng_xgrate), str(porc)]) + '\n')

        
if __name__ == '__main__':
    with open('OINC3.mDBF.subsampled.filtered.convs.pkl', 'rb') as infh:
        convs = pickle.load(infh)
    with open('read2gene.pkl', 'rb') as infh:
        reads2gene = pickle.load(infh)
    numreadspergene, convsPerGene = getPerGene(convs, reads2gene)
    writeConvsPerGene(numreadspergene, convsPerGene, 'test.txt')


    