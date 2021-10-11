#Bioinformatic Analysis of the Conversion Of Nucleotides (BACON)

#Identify genes whose PORC values are different across conditions
#Two possibilities: GLM (readconditions, makePORCdf, getLMEp)
#or
#subsample reads, calculate many PORC values, using Hotelling T2 to identify genes with different PORC distributions across conditions 

import pandas as pd
import sys
import numpy as np
from collections import OrderedDict
import math
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from scipy.stats.distributions import chi2
import warnings
import time


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
            
def getLMEp(sampconddf, porcdf, conditionA, conditionB):
    #Calculate pvalues for genes based on PORC values across conditions using LME model 
    #Delta porc values are calculated as conditionB - condition A

    deltaporcdict = OrderedDict() #{genename : deltaporc} ordered so it's easy to match it up with porcdf
    pvaluedict = OrderedDict() #{genename : pvalue} ordered so it's easy to match it up with q values

    #Get column names in porcdf that are associated with each condition
    conditionAsamps = sampconddf.loc[sampconddf['condition'] == conditionA]
    conditionAsamps = conditionAsamps['sample'].tolist()
    conditionAcolumns = ['porc_' + samp for samp in conditionAsamps]

    conditionBsamps = sampconddf.loc[sampconddf['condition'] == conditionB]
    conditionBsamps = conditionBsamps['sample'].tolist()
    conditionBcolumns = ['porc_' + samp for samp in conditionBsamps]

    print('Condition A samples: ' + (', ').join(conditionAsamps))
    print('Condition B samples: ' + (', ').join(conditionBsamps))

    #Store relationships of conditions and the samples in that condition
    #It's important that this dictionary be ordered because we are going to be iterating through it
    samp_conds = OrderedDict({'condA' : conditionAcolumns, 'condB' : conditionBcolumns})

    #Get a list of all samples
    samps = []
    for cond in samp_conds:
        samps += samp_conds[cond]

    #Iterate through rows, making a dictionary from every row, turning it into a dataframe, then calculating p value
    genecounter = 0
    for index, row in porcdf.iterrows():
        genecounter +=1
        if genecounter % 1000 == 0:
            print('Calculating pvalue for gene {0}...'.format(genecounter))
        
        d = {}
        d['Gene'] = [row['Gene']] * len(samps)
        d['variable'] = samps
        
        values = [] #porc values
        for cond in samp_conds:
            for sample in samp_conds[cond]:
                value = row[sample]
                values.append(value)
        d['value'] = values

        #If there is an NA psi value, we are not going to calculate a pvalue for this gene
        p = None
        if True in np.isnan(values):
            p = np.nan

        conds = []
        for cond in samp_conds:
            conds += [cond] * len(samp_conds[cond])
        condAs = []
        condBs = []
        for cond in conds:
            if cond == 'condA':
                condAs.append(1)
                condBs.append(0)
            elif cond == 'condB':
                condAs.append(0)
                condBs.append(1)
        d['condA'] = condAs #e.g. [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        d['condB'] = condBs #e.g. [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]

        d['samples'] = [x + 1 for x in range(len(samps))]

        #Turn this dictionary into a DataFrame
        rowdf = pd.DataFrame.from_dict(d)

        #delta psi is difference between mean psi of two conditions (cond2 - cond1)
        condAmeanporc = float(format(np.mean(rowdf.query('condA == 1').value.dropna()), '.3f'))
        condBmeanporc = float(format(np.mean(rowdf.query('condB == 1').value.dropna()), '.3f'))
        deltaporc = condBmeanporc - condAmeanporc
        deltaporc = float(format(deltaporc, '.3f'))
        deltaporcdict[row['Gene']] = deltaporc

    #Get LME pvalue, but only if we haven't already determined that the pvalue is NA because we are missing one or more psi values
    #Lots of warnings about convergence, etc. Suppress them.
        if not p:
            with warnings.catch_warnings():
                    warnings.filterwarnings('ignore')

                    #So apparently, some combinations of psi values will give nan p values due to a LinAlgError that arises from a singular
                    #hessian matrix during the fit of the model.  However, the equivalent code in R (nlme::lme) never gives this error, even with
                    #the same data. It's not clear from just looking at the psi values why this is.  However, I found that by varying the
                    #start_params in the fit, this can be avoided. If this is done, the resulting p value always matches what is given in R.
                    #Further, the p value is the same regardless of the start_param.
                    #But it's not clear to me why changing the start_param matters, or what the default is here or with nlme.
                    #So let's try a few starting paramters.  Regardless, this seems to affect a small number of genes (<1%), and it is causing
                    #false negatives because genes that should get p values (may or may not be sig) are getting NA.
                    possible_start_params = [0, 0, 1, -1, 2, -2]
                    numberoftries = -1
                    for param in possible_start_params:
                        #if we already have a pvalue, don't try again
                        if p != None and not np.isnan(p):
                            break
                        #First time through, numberoftries = 0, and we are just using a placeholder startparam (0) here because we aren't even using it.
                        #Gonna use whatever the default is
                        numberoftries += 1
                        try:
                            #actual model
                            md = smf.mixedlm('value ~ condA', data=rowdf, groups='samples', missing='drop')
                            if numberoftries == 0:
                                # REML needs to be false in order to use log-likelihood for pvalue calculation
                                mdf = md.fit(reml=False)
                            elif numberoftries > 0:
                                mdf = md.fit(reml=False, start_params=[param])

                            #null model
                            nullmd = smf.mixedlm('value ~ 1', data=rowdf, groups='samples', missing='drop')
                            if numberoftries == 0:
                                nullmdf = nullmd.fit(reml=False)
                            elif numberoftries > 0:
                                nullmdf = nullmd.fit(reml=False, start_params=[param])

                            #Likelihood ratio
                            LR = 2 * (mdf.llf - nullmdf.llf)
                            p = chi2.sf(LR, df=1)

                        #These exceptions are needed to catch cases where either all psi values are nan (valueerror) or all psi values for one condition are nan (linalgerror)
                        except (ValueError, np.linalg.LinAlgError):
                            p = np.nan

        pvaluedict[row['Gene']] = float('{:.2e}'.format(p))

    #Correct pvalues using BH method, but only using pvalues that are not NA
    pvalues = list(pvaluedict.values())
    pvaluesformultitest = [pval for pval in pvalues if str(pval) != 'nan']
    fdrs = list(multipletests(pvaluesformultitest, method='fdr_bh')[1])
    fdrs = [float('{:.2e}'.format(fdr)) for fdr in fdrs]

    #Now we need to incorporate the places where p = NA into the list of FDRs (also as NA)
    fdrswithnas = []
    fdrindex = 0
    for pvalue in pvalues:
        if str(pvalue) != 'nan':
            fdrswithnas.append(fdrs[fdrindex])
            fdrindex += 1
        elif str(pvalue) == 'nan':
            fdrswithnas.append(np.nan)

    #Add deltaporcs, pvalues, and FDRs to df
    deltaporcs = list(deltaporcdict.values())
    porcdf = porcdf.assign(deltaporc=deltaporcs)
    porcdf = porcdf.assign(pval=pvalues)
    porcdf = porcdf.assign(FDR=fdrswithnas)

    #Write df
    fn = 'porc.pval.txt'
    porcdf.to_csv(fn, sep='\t', index=False, float_format='%.3g', na_rep='NA')

    
    

if __name__ == '__main__':
    sampconddf = readconditions(sys.argv[1])
    #print(sampconddf)
    porcdf = makePORCdf(sys.argv[1], sys.argv[2])
    getLMEp(sampconddf, porcdf, sys.argv[3], sys.argv[4])
