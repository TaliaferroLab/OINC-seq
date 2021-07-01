import pandas as pd
import sys
from functools import reduce
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats.distributions import chi2
import numpy as np
from statsmodels.stats.multitest import multipletests
from collections import OrderedDict
from statsmodels.tools.sm_exceptions import PerfectSeparationError as PSE


def combinesamples(slamdunkouts, samplenames):
    #Given a slew of slam dunk outputs, combine them together into one df
    #slamdunkouts is a comma-separated list of filepaths
    #samplenames is a comma-separated list of names, in the same order as slamdunkouts
    #conditions is a comma-separated list of conditionIDs (max 2), in the same order of slamdunkouts

    dfs = []
    samplenames = samplenames.split(',')
    #conditions = conditions.split(',')
    for idx, sd in enumerate(slamdunkouts.split(',')):
        #skip first 2 header lines
        df = pd.read_csv(sd, sep = '\t', skiprows = 2, header = 0)
        columns = list(df.columns)
        columnstokeep = ['Name', 'G_G', 'G_C', 'G_T']
        columnstodrop = [c for c in columns if c not in columnstokeep]
        df = df.drop(labels = columnstodrop, axis = 1)
        #Combine G_C and G_T
        df = df.assign(G_mut = df['G_C'] + df['G_T'])
        #totalGs
        df = df.assign(totalG = df['G_G'] + df['G_mut'])
        
        #rename columns in preparation for merging
        columns = list(df.columns)
        untouchablecolumnnames = ['Name']
        samplename = samplenames[idx]
        columns = [c + ';' + samplename if c not in untouchablecolumnnames else c for c in columns]
        df.columns = columns
        df = df.drop_duplicates(ignore_index = True) #somehow duplicate transcripts?
        dfs.append(df)

    bigdf = reduce(lambda x, y: pd.merge(x, y, on = ['Name']), dfs)
    bigdf = bigdf.drop_duplicates(ignore_index = True)
    #Remove any row with NA value
    bigdf = bigdf.dropna(axis = 0, how = 'any')
    
    return bigdf

def classifysamples(samplenames, conditions):
    #samplenames is a comma-separated list of names
    #conditions is a comma-separated list of conditionIDs (max 2), in the same order of samplenames
    samplenames = samplenames.split(',')
    conditions = conditions.split(',')
    d = dict(zip(samplenames, conditions))
    
    return d

def doglm(bigdf, sampconds):
    genecounter = 0
    ps = OrderedDict() #{genename : p}
    mintotalGs = OrderedDict() #{genename : min number of total Gs across all samples}
    meantotalGs = OrderedDict() #{genename : mean number of total Gs across all samples}
    condArates = OrderedDict() #{genename : condA mean mutation rate}
    condBrates = OrderedDict() #{genename : condA mean mutation rate}
    ratelog2fc = OrderedDict() #{genename : log2fc in mutation rates (B/A)}
    #sampconds is dict of {samplename : condition}
    for index, row in bigdf.iterrows():
        genecounter +=1
        if genecounter % 1000 == 0:
            print('Gene {0}...'.format(genecounter))
        samples = list(sampconds.keys())
        G_Gs = []
        G_Cs = []
        G_Ts = []
        G_muts = []
        totalGs = []
        conds = []
        for sample in samples:
            G_G = row['G_G;{0}'.format(sample)]
            G_C = row['G_C;{0}'.format(sample)]
            G_T = row['G_T;{0}'.format(sample)]
            G_mut = row['G_mut;{0}'.format(sample)]
            totalG = row['totalG;{0}'.format(sample)]
            cond = sampconds[sample]
            G_Gs.append(G_G)
            G_Cs.append(G_C)
            G_Ts.append(G_T)
            G_muts.append(G_mut)
            totalGs.append(totalG)
            conds.append(cond)

        mintotalGs[row['Name']] = min(totalGs)
        meantotalGs[row['Name']] = np.mean(totalGs)
        d = {'G_G' : G_Gs, 'G_C' : G_Cs, 'G_T' : G_Ts, 'G_mut' : G_muts, 'totalG' : totalGs, 'cond' : conds}
        rowdf = pd.DataFrame.from_dict(d)
        totalGs = rowdf['totalG'].tolist()
        totalmuts = rowdf['G_mut'].tolist()
        minG = min(totalGs)
        minmut = min(totalmuts)
        #Implement totalG count filter
        if minG < 100:
            p = np.nan

        else:           
            try:
                #GLM
                mod_real = smf.glm('G_mut + G_G ~ cond', family = sm.families.Binomial(), data = rowdf).fit()
                mod_null = smf.glm('G_mut + G_G ~ 1', family = sm.families.Binomial(), data = rowdf).fit()

                #Likelihood ratio test
                logratio = (mod_real.llf - mod_null.llf) * 2
                p = round(chi2.sf(logratio, df = 1), 4)
            #If all mutation counts in one condition are 0, this causes a problem called Perfect Separation Error
            #Interestingly, this is not triggered if there is only 1 replicate in a condition
            except PSE:
                p = np.nan


        ps[row['Name']] = p

        #Calculate mean conversion rates for each condition and the difference between them
        conds = sorted(list(set(conds))) #conds are alphabetically sorted. condA is the first one, condB is the second one
        for idx, cond in enumerate(conds):
            conddf = rowdf[rowdf['cond'] == cond]
            rates = []
            for condi, condr in conddf.iterrows():
                try:
                    rate = condr['G_mut'] / condr['totalG']
                except ZeroDivisionError:
                    rate = np.nan

                rates.append(rate)

            meanrate = np.mean(rates)
            if idx == 0:
                condArates[row['Name']] = meanrate
                condArate = meanrate
            elif idx == 1:
                condBrates[row['Name']] = meanrate
                condBrate = meanrate

        pc = 1e-6
        log2fc = np.log2((condBrate + pc) / (condArate + pc))
        ratelog2fc[row['Name']] = log2fc

    #Correct pvalues using BH method, but only using pvalues that are not NA
    pvalues = list(ps.values())
    pvaluesformultitest = [pval for pval in pvalues if str(pval) != 'nan']
    fdrs = list(multipletests(pvaluesformultitest, method = 'fdr_bh')[1])
    fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

    #Now we need to incorporate the places where p = NA into the list of FDRs (also as NA)
    fdrswithnas = []
    fdrindex = 0
    for pvalue in pvalues:
        if str(pvalue) != 'nan':
            fdrswithnas.append(fdrs[fdrindex])
            fdrindex +=1
        elif str(pvalue) == 'nan':
            fdrswithnas.append(np.nan)

    genes = bigdf['Name'].tolist()
    outd = {'Gene' : genes, 'minGcount' : list(mintotalGs.values()), 'meanGcount' : list(meantotalGs.values()), '{0}mutrate'.format(conds[0]) : list(condArates.values()), '{0}mutrate'.format(conds[1]) : list(condBrates.values()), 'log2fc' : list(ratelog2fc.values()), 'p' : list(ps.values()), 'FDR' : fdrswithnas}
    outdf = pd.DataFrame.from_dict(outd)

    #format output columns
    formats = {'minGcount': '{:d}', 'meanGcount': '{:.2f}', '{0}mutrate'.format(conds[0]): '{:.2e}', '{0}mutrate'.format(conds[1]): '{:.2e}', 'log2fc': '{:.2f}', 'p': '{:.3f}', 'FDR': '{:.3f}'}
    for col, f in formats.items():
        outdf[col] = outdf[col].map(lambda x : f.format(x))
    
    outdf.to_csv('glm.txt', sep = '\t', header = True, index = False, na_rep = 'NA')      


bigdf = combinesamples(sys.argv[1], sys.argv[2])
sampconds = classifysamples(sys.argv[2], sys.argv[3])
doglm(bigdf, sampconds)
