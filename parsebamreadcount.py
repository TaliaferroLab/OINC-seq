import sys
import os
import argparse

#Take the output of bam-readcount and parse to get a table of nt fractions at each position

def parsebrc(bamreadcountout, covthresh, outfile):
    d = {} #{chrm : {position : [ref, depth, Afrac, Tfrac, Gfrac, Cfrac]}}

    with open(bamreadcountout, 'r') as infh:
        for line in infh:
            line = line.strip().split('\t')
            chrm = line[0]
            position = int(line[1])
            ref = line[2]
            depth = int(line[3])
            if depth < covthresh:
                continue

            if chrm not in d:
                d[chrm] = {}
            d[chrm][position] = [ref, depth]
            ntfracs = {} #{nt : frac of reads}
            for f in line[4:]:
                fsplit = f.split(':')
                nt = fsplit[0]
                if nt not in ['A', 'C', 'T', 'G']:
                    continue
                ntcount = int(fsplit[1])
                ntfrac = ntcount / depth
                ntfrac = f'{ntfrac:.2e}'
                ntfrac = float(ntfrac)
                if nt == ref:
                    ntfrac = 'NA'
                ntfracs[nt] = ntfrac

            d[chrm][position].extend([ntfracs['A'], ntfracs['T'], ntfracs['G'], ntfracs['C']])

    with open(outfile, 'w') as outfh:
        outfh.write(('\t').join(['chrm', 'position', 'ref', 'depth', 'Afrac', 'Tfrac', 'Gfrac', 'Cfrac']) + '\n')
        for chrm in d:
            for position in d[chrm]:
                ref, depth, afrac, tfrac, gfrac, cfrac = d[chrm][position]
                outfh.write(('\t').join([chrm, str(position), ref, str(depth), str(afrac), str(tfrac), str(gfrac), str(cfrac)]) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bamreadcountout', type = str, help = 'Output of bam-readcount.')
    parser.add_argument('--mindepth', type = int, help = 'Minimum depth for a position to be considered.')
    parser.add_argument('--outfile', type = str, help = 'Output file.')
    args = parser.parse_args()

    parsebrc(args.bamreadcountout, args.mindepth, args.outfile)

    