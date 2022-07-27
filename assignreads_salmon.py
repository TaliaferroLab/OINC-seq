import os
import sys
import pysam
import pickle
import gffutils
import numpy as np


#Take in a dictionary of {readid : conversions} (made by getmismatches.py) and a postmaster-enhanced bam (made by alignAndQuant.py).
#First, construct dictionary of {readid : {txid : fractional assignment}}. Then, combining this dictionary with the previous one,
#count the number of conversions associated with each transcript. Finally (and I guess optionally), using a genome annotation file,
#collapse transcript level conversion counts to gene-level conversion counts.

def getpostmasterassignments(postmasterbam):
    #Given a postmaster-produced bam, make a dictionary of the form {readid : {txid : fractional assignment}}
    #It looks like in a postmaster bam that paired end reads are right after each other and are always
    #given the same fractional assignments. This means we can probably just consider R1 reads.

    pprobs = {} #{readid : {txid : pprob}}

    with pysam.AlignmentFile(postmasterbam, 'r') as bamfh:
        for read in bamfh.fetch(until_eof = True):
            if read.is_read2:
                continue
            readid = read.query_name
            tx = read.reference_name
            pprob = read.get_tag(tag='ZW')
            if readid not in pprobs:
                pprobs[readid] = {}
            pprobs[readid][tx] = pprob

    return pprobs

def assigntotxs(pprobs, convs):
    #Intersect posterior probabilities of read assignments to transcripts with conversion counts of those reads.
    #The counts assigned to a tx by a read are scaled by the posterior probability that a read came from that transcript.

    #pprobs = #{readid : {txid : pprob}}
    #produced from getpostmasterassignments()
    #convs = #{readid : {a_a : 200, a_t : 1, etc.}}
    print('Finding transcript assignments for {0} reads.'.format(len(convs)))
    readswithoutassignment = 0 #number of reads which exist in convs but not in pprobs (i.e. weren't assigned to a transcript by salmon)

    txconvs = {}  # {txid : {a_a : 200, a_t : 1, etc.}}

    for readid in pprobs:

        try:
            readconvs = convs[readid]
        except KeyError: #we couldn't find this read in convs
            readswithoutassignment +=1
            continue

        for txid in pprobs[readid]:
            if txid not in txconvs:
                txconvs[txid] = {}
            pprob = pprobs[readid][txid]
            for conv in readconvs:
                scaledconv = readconvs[conv] * pprob
                if conv not in txconvs[txid]:
                    txconvs[txid][conv] = scaledconv
                else:
                    txconvs[txid][conv] += scaledconv

    readswithtxs = len(convs) - readswithoutassignment
    pct = round(readswithtxs / len(convs), 2) * 100
    print('Found transcripts for {0} of {1} reads ({2}%).'.format(readswithtxs, len(convs), pct))

    return txconvs

def collapsetogene(txconvs, gff):
    #Collapse tx-level count measurements to gene level.
    #Need to relate transcripts and genes. Do that with the supplied gff annotation.
    #txconvs = {txid : {a_a : 200, a_t : 1, etc.}}
    
    tx2gene = {} #{txid : geneid}
    geneid2genename = {} #{geneid : genename}
    geneconvs = {}  # {geneid : {a_a : 200, a_t : 1, etc.}}

    print('Indexing gff..')
    gff_fn = gff
    db_fn = os.path.abspath(gff_fn) + '.db'
    if os.path.isfile(db_fn) == False:
        gffutils.create_db(gff_fn, db_fn, merge_strategy='merge', verbose=True)
    print('Done indexing!')

    db = gffutils.FeatureDB(db_fn)
    genes = db.features_of_type('gene')

    print('Connecting transcripts and genes...')
    for gene in genes:
        geneid = str(gene.id).split('.')[0] #remove version numbers
        genename = gene.attributes['gene_name'][0]
        geneid2genename[geneid] = genename
        for tx in db.children(gene, featuretype = 'transcript'):
            txid = str(tx.id).split('.')[0]
            tx2gene[txid] = geneid
    print('Done!')

    allgenes = list(set(tx2gene.values()))

    #Initialize geneconvs dictionary
    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n']

    for gene in allgenes:
        geneconvs[gene] = {}
        for conv in possibleconvs:
            geneconvs[gene][conv] = 0

    for tx in txconvs:
        try:
            gene = tx2gene[tx]
        except KeyError:
            print('WARNING: transcript {0} doesn\'t belong to a gene in the supplied annotation.'.format(tx))
            continue
        convs = txconvs[tx]
        for conv in convs:
            convcount = txconvs[tx][conv]
            geneconvs[gene][conv] += convcount

    return tx2gene, geneid2genename, geneconvs

def readspergene(quantsf, tx2gene):
    #Get the number of reads assigned to each tx. This can simply be read from the salmon quant.sf file.
    #Then, sum read counts across all transcripts within a gene. 
    #Transcript and gene relationships were derived by collapsetogene().

    txcounts = {} #{txid : readcounts}
    genecounts = {} #{geneid : readcounts}

    with open(quantsf, 'r') as infh:
        for line in infh:
            line = line.strip().split('\t')
            if line[0] == 'Name':
                continue
            txid = line[0].split('.')[0] #remove tx id version in the salmon quant.sf if it exists
            counts = float(line[4])
            txcounts[txid] = counts

    allgenes = list(set(tx2gene.values()))
    for gene in allgenes:
        genecounts[gene] = 0

    for txid in txcounts:
        try:
            geneid = tx2gene[txid]
        except KeyError: #maybe the salmon tx id have version numbers
            txid = txid.split('.')[0]
            geneid = tx2gene[txid]
        
        genecounts[geneid] += txcounts[txid]

    return genecounts


def writeOutput(sampleparams, geneconvs, genecounts, geneid2genename, outfile, use_g_t, use_g_c):
    #Write number of conversions and readcounts for genes.
    possibleconvs = [
        'a_a', 'a_t', 'a_c', 'a_g', 'a_n',
        'g_a', 'g_t', 'g_c', 'g_g', 'g_n',
        'c_a', 'c_t', 'c_c', 'c_g', 'c_n',
        't_a', 't_t', 't_c', 't_g', 't_n']

    with open(outfile, 'w') as outfh:
        #Write arguments for this pigpen run
        for arg in sampleparams:
            outfh.write('#' + arg + '\t' + str(sampleparams[arg]) + '\n')
        #total G is number of ref Gs encountered
        #convG is g_t + g_c (the ones we are interested in)
        outfh.write(('\t').join(['GeneID', 'GeneName', 'numreads'] + possibleconvs + [
                    'totalG', 'convG', 'convGrate', 'G_Trate', 'G_Crate', 'porc']) + '\n')
        genes = sorted(geneconvs.keys())

        for gene in genes:
            genename = geneid2genename[gene]
            numreads = genecounts[gene]
            convcounts = []
            c = geneconvs[gene]
            for conv in possibleconvs:
                convcount = c[conv]
                convcounts.append(convcount)

            convcounts = ['{:.2f}'.format(x) for x in convcounts]

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
            if type(numreads) == float:
                numreads = '{:.2f}'.format(numreads)
            if type(convG) == float:
                convG = '{:.2f}'.format(convG)
            if type(totalG) == float:
                totalG = '{:.2f}'.format(totalG)
            if type(convGrate) == float:
                convGrate = '{:.2e}'.format(convGrate)
            if type(g_trate) == float:
                g_trate = '{:.2e}'.format(g_trate)
            if type(g_crate) == float:
                g_crate = '{:.2e}'.format(g_crate)
            if type(porc) == np.float64:
                porc = '{:.3f}'.format(porc)

            outfh.write(('\t').join([gene, genename, str(numreads)] + convcounts + [str(totalG), str(convG), str(convGrate), str(g_trate), str(g_crate), str(porc)]) + '\n')

    

    

if __name__ == '__main__':
    print('Getting posterior probabilities from salmon alignment file...')
    pprobs = getpostmasterassignments(sys.argv[1])
    print('Done!')
    print('Loading conversions from pickle file...')
    with open(sys.argv[2], 'rb') as infh:
        convs = pickle.load(infh)
    print('Done!')
    print('Assinging conversions to transcripts...')
    txconvs = assigntotxs(pprobs, convs)
    print('Done!')

    tx2gene, geneid2genename, geneconvs = collapsetogene(txconvs, sys.argv[3])
    genecounts = readspergene(sys.argv[4], tx2gene)
    writeOutput(geneconvs, genecounts, geneid2genename, sys.argv[5], True, True)