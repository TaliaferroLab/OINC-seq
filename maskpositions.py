#Given a bed of positions to mask, extract them and store them in a dictionary
#so they can be given to getmismatches.py.
#Beds are 0-based half-open, and we want to give getmismatches 0-based coordinates for easy
#interfacing with pysam.

def readmaskbed(maskbed):
    maskdict = {} #{chrm : (set of positions to mask)}
    with open(maskbed, 'r') as infh:
        for line in infh:
            line = line.strip().split('\t')
            chrm = line[0]
            start = int(line[1])
            end = int(line[2])
            interval = list(range(start, end))
            if chrm not in maskdict:
                maskdict[chrm] = []
            maskdict[chrm].extend(interval)


    #Remove duplicates
    for chrm in maskdict:
        coords = set(maskdict[chrm])
        maskdict[chrm] = coords

    #Tell us how many positions we are masking
    totalmask = 0
    for chrm in maskdict:
        totalmask += len(maskdict[chrm])

    print('Manually masking {0} positions.'.format(totalmask))

    return maskdict