# OINC-seq <br/> <br/>Detecting oxidative marks on RNA using high-throughput sequencing

## Overview

OINC-seq (Oxidation-Induced Nucleotide Conversion sequencing) is a sequencing technology that allows the direction of oxidative marks on RNA molecules. Because guanosine has the lowest redox potential of any of the ribonucleosides, it is the one most likely to be affected by oxidation. When this occurs, guanosine is turned into 8-oxoguanosine (8-OG). A previous [study](https://pubs.acs.org/doi/10.1021/acs.biochem.7b00730) found that when reverse transcriptase encounters guanosine oxidation products, it can misinterpret 8-OG as either T or C. Therefore, to detect these oxidative marks, one can look for G -> T and G -> C conversions in RNAseq data.

To detect and quantify these conversions, we have created software called **PIGPEN** (Pipeline for Identification of Guanosine Positions Erroneously Notated).

PIGPEN takes in alignment files (bam), ideally made with [STAR](https://github.com/alexdobin/STAR). Single and paired-end reads are supported, although paired-end reads are preferred (for reasons that will become clear later). To minimize the contribution of positions that appear as mutations due to non-ideal alignments, PIGPEN only considers uniquely aligned reads (mapping quality == 255). For now, it is required that paired-end reads be stranded, and that read 1 correspond to the sense strand. This is true for most, but not all, modern RNAseq library preparation protocols.

                        ,-,-----,
        PIGPEN     **** \ \ ),)`-'
                  <`--'> \ \`
                  /. . `-----,
        OINC! >  ('')  ,      @~
                  `-._,  ___  /
    -|-|-|-|-|-|-|-| (( / (( / -|-|-|
    |-|-|-|-|-|-|-|- ''' ''' -|-|-|-
    -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

    Pipeline for Identification
       Of Guanosine Positions
        Erroneously Notated

## Requirements

PIGPEN has the following prerequisites:

- python >= 3.6
- samtools >= 1.13
- varscan >= 2.4.4
- bcftools >= 1.13
- pybedtools >= 0.8.2
- pysam >= 0.16
- numpy >= 1.21
- pandas >= 1.3.3
- bamtools >= 2.5.1
- bedtools >= 2.30.0

BACON has the following prerequisites:

- python >= 3.6
- statsmodels >= 0.13.2
- numpy >= 1.21
- rpy2 >= 3.4.5
- R >= 4.1

## Installation

For now, installation can be done by cloning this repository. As PIPGEN matures, we will work towards getting this package on [bioconda](https://bioconda.github.io/).

## SNPs

8-OG-induced conversions are rare, and this rarity makes it imperative that contributions from conversions that are not due to oxidation are minimized. A major source of apparent conversions is SNPs. It is therefore advantageous to find and mask SNPs in the data.

PIGPEN performs this by using [varscan](http://varscan.sourceforge.net/using-varscan.html) to find SNP positions. These locations are then excluded from all future analyses. Varscan parameters are controled by the PIGPEN parameters `--SNPcoverage` and `--SNPfreq` that control the depth and frequency required to call a SNP. We recommend being aggressive with these parameters. We often set them to 20 and 0.02, respectively.

PIGPEN performs this SNP calling on control alignment files (`--controlBams`) in which the intended oxidation did not occur. PIGPEN will use the union of all SNPs found in these files for masking. Whether or not to call SNPs at all (you probably should) is controlled by `--useSNPs`.

This process can be time consuming. At the end, a file called **merged.vcf** is created in the current working directory. If this file is present, PIGPEN will assume that it should be used for SNP masking, allowing the process of identifying SNPs to be skipped.

## Filtering alignments

Because the process of finding nucleotide conversions can take a long time, PIGPEN first filters the reads, keeping only those that overlap with any feature in a supplied bed file (`--geneBed`). This process can use multiple processors (`--nproc`) to speed it up, and requires another file (`--chromsizes`). This file is a 2 column, tab-delimited text file where column 1 is the reference (chromosome) names for the references present in the alignment file, and column 2 is the integer size of that reference. If a fasta file and fasta index for the genome exists, this file can be made using `cut -f 1,2 genome.fa.fai`.

## Quantifying conversions

PIGPEN then identifies conversions in reads. This can be done using multiple processors (`--nproc`). In order to minimize the effect of sequencing error, PIGPEN only considers positions for which the sequencing quality was at least 30. There are two important flags to consider here.

First, `--onlyConsiderOverlap` requires that the same conversion be observed in both reads of a mate pair. Positions interrogated by only one read are not considered. This can improve accuracy. True oxidation-induced conversions are rare. Rare enough that sequencing errors can cause a problem. Requiring that a conversion be present in both reads minimizes the effect of sequencing errors. If the fragment sizes for a library are especially large relative to the read length, the number of positions interrogated by both mates will be small.

Second, `--requireMultipleConv` requires that there be at least two G -> C / G -> T conversions in a read pair in order for those conversions to be recorded. The rationale here is again to reduce the contribution of background, non-oxidation-related conversions. Background conversions should be distributed relatively randomly across reads. However, due to the spatial nature of the oxidation reaction, oxidation-induced conversions should be more clustered into specific reads. Therefore, requiring at least two conversions can increase specificity. In practice, this works well if the data is very deep or concentrated on a small number of targets. When dealing with transcriptome-scale data, this flag often reduces the number of observed conversions to an unacceptably low level.

## Assigning reads to genes

For now, PIGPEN using `bedtools` and a supplied bed file of gene locations (`--geneBed`) to assign individual reads to genes. We are working on improvements in this area.

## Calculating the number of conversions per gene

After identifying the conversions present in each read and the cognate gene for each read, the number of conversions for each gene is calculated. We have observed that the overall rate of conversions (not just G -> T + G -> C, but all conversions) can vary signficantly from sample to sample, presumably due to a technical effect in library preparation. For this reason, PIGPEN calculates **PORC** (Proportion of Relevant Conversions) values. This is the log2 ratio of the relevant conversion rate ([G -> T + G -> C] / total number of reference G encountered) to the overall conversion rate (total number of all conversions / total number of positions interrogated). PORC therefore normalizes to the overall rate of conversions, removing this technical effect.

PIGPEN can use G -> T conversions, G -> C conversions, or both when calculating PORC values. This behavior is controlled by supplying the options `--use_g_t` and `--use_g_c`.

## Statistical framework for comparing gene-level PORC values across conditions

We could simply compare PORC values across conditions, but with that approach we lose information about the number of counts (conversions) that went into the PORC calculation.

For each gene, PIGPEN calculates the number of relevant conversions (G -> T + G -> C) as well as all other conversions encountered. Each gene therefore ends up with a 2x2 contingency table of the following form:

|| converted | not converted |
| ----------------|---------------|-------- |
| G -> C or G -> T | a | b | 
| other conversion | c | d | 

We then want to compare groups (replicates) of contingency tables across conditions. BACON (Bioinformatic Analysis of the Conversion of Nucleotides) performs this comparison using a binomial linear mixed-effects model. Replicates are modeled as random effects.

`full model = conversions ~ nucleotide + condition + nucleotide:condition + (1 | replicate)`

`null model = conversions ~ nucleotide + condition + (1 | replicate)`

The two models are then compared using a likelihood ratio test.

As input, BACON takes a tab-delimited, headered file of the following form with one row per sample:

| file | sample | condition |
| -----|--------|-----------|
| /path/to/pigpen/output | sample name | condition ID|
