#Bioinformatic Analysis of the Conversion Of Nucleotides (BACON)

#Use a linear model to identify genes whose rate of G-conversions changes across conditions
#Going to end up with a lot of contingency tables:

#       converted   notConverted
#G
#nonG

#G conversions are only allowed to be G->T and G->C
#Compare groups of contingency tables across conditions
import pandas as pd
import sys
import numpy as np
import itertools
from statsmodels.stats.multitest import multipletests
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, Formula, FloatVector
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.rinterface import RRuntimeWarning
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
import warnings

#Need r-base, r-stats, r-lme4

#supress RRuntimeWarnings
warnings.filterwarnings('ignore', category = RRuntimeWarning)
rpy2_logger.setLevel(logging.ERROR) #suppresses R messages to console



def readconditions(samp_conds_file):
    #Read in a three column file of oincoutputfile / sample / condition
    #Return a pandas df 
    sampconddf = pd.read_csv(samp_conds_file, sep = '\t', index_col = False, header = 0)

    return sampconddf

def makePORCdf(samp_conds_file, minreads):
    #Make a dataframe of PORC values for all samples
    #start with GENE...SAMPLE...READCOUNT...PORC
    #then make a wide version that is GENE...SAMPLE1READCOUNT...SAMPLE1PORC...SAMPLE2READCOUNT...SAMPLE2PORC

    #Only keep genes that are present in EVERY pigpen file

    minreads = int(minreads)
    dfs = [] #list of dfs from individual pigpen outputs
    genesinall = set() #genes present in every pigpen file
    with open(samp_conds_file, 'r') as infh:
        for line in infh:
            line = line.strip().split('\t')
            if line[0] == 'file':
                continue
            pigpenfile = line[0]
            sample = line[1]
            df = pd.read_csv(pigpenfile, sep = '\t', index_col = False, header = 0)
            dfgenes = df['Gene'].tolist()
            samplecolumn = [sample] * len(dfgenes)
            df = df.assign(sample = samplecolumn)

            if not genesinall: #if there are no genes in there (this is the first file)
                genesinall = set(dfgenes)
            else:
                genesinall = genesinall.intersection(set(dfgenes))
            
            columnstokeep = ['Gene', 'sample', 'numreads', 'porc'] 
            df = df[columnstokeep]
            dfs.append(df)

    #for each df, filter keeping only the genes present in every df (genesinall)
    #Somehow there are some genes whose name in NA
    if np.nan in genesinall:
        genesinall.remove(np.nan)
    dfs = [df.loc[df['Gene'].isin(genesinall)] for df in dfs]

    #concatenate (rbind) dfs together
    df = pd.concat(dfs)
    
    #turn from long into wide
    df = df.pivot_table(index = 'Gene', columns = 'sample', values = ['numreads', 'porc']).reset_index()
    #flatten multiindex column names
    df.columns = ["_".join(a) if '' not in a else a[0] for a in df.columns.to_flat_index()]
    
    #Filter for genes with at least minreads in every sample
    #get columns with numreads info
    numreadsColumns = [col for col in df.columns if 'numreads' in col]
    #Get minimum in those columns
    df['minreadcount'] = df[numreadsColumns].min(axis = 1)
    #Filter for rows with minreadcount >= minreads
    print('Filtering for genes with at least {0} reads in every sample.'.format(minreads))
    df = df.loc[df['minreadcount'] >= minreads]
    print('{0} genes have at least {1} reads in every sample.'.format(len(df), minreads))
    #We also don't want rows with inf/-inf PORC values
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(how= 'any')
    #Return a dataframe of just genes and PORC values
    columnstokeep = ['Gene'] + [col for col in df.columns if 'porc' in col]
    df = df[columnstokeep]

    return df

def calcDeltaPORC(porcdf, sampconds, conditionA, conditionB):
    #Given a porc df from makePORCdf, add deltaporc values.

    deltaporcs = []
    sampconddf = readconditions(sampconds)

    #Get column names in porcdf that are associated with each condition
    conditionAsamps = sampconddf.loc[sampconddf['condition'] == conditionA]
    conditionAsamps = conditionAsamps['sample'].tolist()
    conditionAcolumns = ['porc_' + samp for samp in conditionAsamps]

    conditionBsamps = sampconddf.loc[sampconddf['condition'] == conditionB]
    conditionBsamps = conditionBsamps['sample'].tolist()
    conditionBcolumns = ['porc_' + samp for samp in conditionBsamps]

    print('Condition A samples: ' + (', ').join(conditionAsamps))
    print('Condition B samples: ' + (', ').join(conditionBsamps))

    for index, row in porcdf.iterrows():
        condAporcs = []
        condBporcs = []
        for col in conditionAcolumns:
            porc = row[col]
            condAporcs.append(porc)
        for col in conditionBcolumns:
            porc = row[col]
            condBporcs.append(porc)

        condAporcs = [x for x in condAporcs if x != np.nan]
        condBporcs = [x for x in condBporcs if x != np.nan]
        condAporc = np.mean(condAporcs)
        condBporc = np.mean(condBporcs)
        deltaporc = condBporc - condAporc
        deltaporc = float(format(deltaporc, '.3f'))
        deltaporcs.append(deltaporc)

    porcdf = porcdf.assign(deltaPORC = deltaporcs)
    
    return porcdf

def makeContingencyTable(row):
    #Given a row from a pigpen df, return a contingency table of the form
    #[[convG, nonconvG], [convnonG, nonconvnonG]]

    convG = row['g_t'] + row['g_c']
    nonconvG = row['g_g']
    convnonG = row['a_t'] + row['a_c'] + row['a_g'] + row['c_a'] + row['c_t'] + row['c_g'] + row['t_a'] + row['t_c'] + row['t_g']
    nonconvnonG = row['c_c'] + row['t_t'] + row['a_a']

    conttable = [[convG, nonconvG], [convnonG, nonconvnonG]]

    return conttable


def calculate_nested_f_statistic(small_model, big_model):
    #Given two fitted GLMs, the larger of which contains the parameter space of the smaller, return the F Stat and P value corresponding to the larger model adding explanatory power
    #No anova test for GLM models in python
    #this is a workaround
    #https://stackoverflow.com/questions/27328623/anova-test-for-glm-in-python
    addtl_params = big_model.df_model - small_model.df_model
    f_stat = (small_model.deviance - big_model.deviance) / \
        (addtl_params * big_model.scale)
    df_numerator = addtl_params
    # use fitted values to obtain n_obs from model object:
    df_denom = (big_model.fittedvalues.shape[0] - big_model.df_model)
    p_value = stats.f.sf(f_stat, df_numerator, df_denom)
    return (f_stat, p_value)

def getgenep(geneconttable):
    #Given a gene-level contingency table of the form:
    #[condAtables, condBtables], where each individual sample is of the form
    #[[convG, nonconvG], [convnonG, nonconvnonG]],
    #run glm either including or excluding condition term
    #using likelihood ratio of the two models and chisq test, return p value

    #Turn gene-level contingency tables into df of form
    #convcount  nonconvcount    condition   nuc sample
    nCondAsamples = len(geneconttable[0])
    nCondBsamples = len(geneconttable[1])

    # e.g. ['condA', 'condA', 'condB', 'condB']
    cond = ['condA'] * (nCondAsamples * 4) + ['condB'] * (nCondBsamples * 4)
    nuc = ['G', 'G', 'nonG', 'nonG'] * nCondAsamples + ['G', 'G', 'nonG',
                                                        'nonG'] * nCondBsamples  # e.g. ['G', 'G', 'nonG', 'nonG', ...]
    conv = ['yes', 'no', 'yes', 'no'] * nCondAsamples + ['yes', 'no',
                                                         'yes', 'no'] * nCondBsamples  # e.g. ['yes', 'no', 'yes', 'no', ...]
    samples = ['sample' + str(x + 1) for x in range(nCondAsamples + nCondBsamples)]
    samples = list(itertools.repeat(samples, 4))
    samples = list(itertools.chain.from_iterable(samples))
    samples = sorted(samples)
    samplenumber = [x + 1 for x in range(len(cond))]

    #to get counts, flatten the very nested list of lists that is geneconttable
    a = list(itertools.chain.from_iterable(geneconttable))
    b = list(itertools.chain.from_iterable(a))
    counts = list(itertools.chain.from_iterable(b))

    d = {'cond': cond, 'nuc': nuc, 'conv': conv,
         'counts': counts, 'sample' : samples}
    df = pd.DataFrame.from_dict(d)
    #Reshape table to get individual columns for converted and nonconverted nts
    df2 = df.pivot_table(index = ['cond', 'nuc', 'sample'], columns = 'conv', values = 'counts').reset_index()

    pandas2ri.activate()

    fmla = 'cbind(yes, no) ~ nuc + cond + nuc:cond + (1 | sample)'
    nullfmla = 'cbind(yes, no) ~ nuc + cond + (1 | sample)'

    fullfit = lme4.glmer(formula=fmla, family=stats.binomial, data=df2)
    reducedfit = lme4.glmer(formula=nullfmla, family=stats.binomial, data=df2)

    logratio = (stats.logLik(fullfit)[0] - stats.logLik(reducedfit)[0]) * 2
    pvalue = stats.pchisq(logratio, df=2, lower_tail=False)[0]
    #format decimal
    pvalue = float('{:.2e}'.format(pvalue))
    
    return pvalue

def multihyp(pvalues):
    #given a dictionary of {gene : pvalue}, perform multiple hypothesis correction

    #remove genes with p value of NA
    cleanedp = {}
    for gene in pvalues:
        if not np.isnan(pvalues[gene]):
            cleanedp[gene] = pvalues[gene]

    cleanedps = list(cleanedp.values())
    fdrs = multipletests(cleanedps, method = 'fdr_bh')[1]
    fdrs = dict(zip(list(cleanedp.keys()), fdrs))

    correctedps = {}
    for gene in pvalues:
        if np.isnan(pvalues[gene]):
            correctedps[gene] = np.nan
        else:
            correctedps[gene] = fdrs[gene]

    return correctedps


def getpvalues(samp_conds_file, conditionA, conditionB):
    #each contingency table will be: [[convG, nonconvG], [convnonG, nonconvnonG]]
    #These will be stored in a dictionary: {gene : [condAtables, condBtables]}
    conttables = {}

    pvalues = {} #{gene : pvalue}

    nsamples = 0
    with open(samp_conds_file, 'r') as infh:
        for line in infh:
            line = line.strip().split('\t')
            if line[0] == 'file':
                continue
            nsamples +=1
            pigpenfile = line[0]
            sample = line[1]
            condition = line[2]
            df = pd.read_csv(pigpenfile, sep = '\t', index_col = False, header=0)
            for idx, row in df.iterrows():
                conttable = makeContingencyTable(row)
                gene = row['Gene']
                if gene not in conttables:
                    conttables[gene] = [[], []]
                if condition == conditionA:
                    conttables[gene][0].append(conttable)
                elif condition == conditionB:
                    conttables[gene][1].append(conttable)

    genecounter = 0
    for gene in conttables:
        genecounter +=1
        if genecounter % 1000 == 0:
            print('Getting p value for gene {0}...'.format(genecounter))
        geneconttable = conttables[gene]
        #Only calculate p values for genes present in every porc file
        gene_porcfiles = len(geneconttable[0]) + len(geneconttable[1])
        if nsamples == gene_porcfiles:
            try:
                p = getgenep(geneconttable)
            except RRuntimeError:
                p = np.nan
        else:
            p = np.nan
        
        pvalues[gene] = p

    correctedps = multihyp(pvalues)
    
    pdf = pd.DataFrame.from_dict(pvalues, orient = 'index', columns = ['pval'])
    fdrdf = pd.DataFrame.from_dict(correctedps, orient = 'index', columns = ['FDR'])
    
    pdf = pd.merge(pdf, fdrdf, left_index = True, right_index = True).reset_index().rename(columns = {'index' : 'Gene'})

    return pdf

def formatporcDF(porcdf):
    #Format floats in all porcDF columns
    formats = {'deltaPORC': '{:.3f}', 'pval': '{:.3e}', 'FDR': '{:.3e}'}
    c = porcdf.columns.tolist()
    c = [x for x in c if 'porc_' in x]
    for x in c:
        formats[x] = '{:.3f}' #all porc_SAMPLE columns
    for col, f in formats.items():
        porcdf[col] = porcdf[col].map(lambda x: f.format(x))

    return porcdf

if __name__ == '__main__':
    utils = importr('utils')
    lme4 = importr('lme4')
    base = importr('base')
    stats = importr('stats')

    #Make df of PORC values
    porcdf = makePORCdf(sys.argv[1], 100)
    #Add delta porc values
    porcdf = calcDeltaPORC(porcdf, sys.argv[1], 'mDBF', 'pDBF')
    #Get p values and corrected p values
    pdf = getpvalues(sys.argv[1], 'mDBF', 'pDBF')
    #Add p values and FDR
    porcdf = pd.merge(porcdf, pdf, on = ['Gene'])
    #Format floats
    porcdf = formatporcDF(porcdf)

    porcdf.to_csv('porc.txt', sep='\t', index = False)
