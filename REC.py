#Input: file,  Its head must contain a line of cancers. It must also contain a list of Mirna to be tested.  then tab deliminiated canger, mirna, linreg results: pval.
#output: rec score or recurrence scores of cancer.  Using a relative rank to make an even distribution the sum of these can approximate a chisquared test.  We use this for pval for the null hypothesis and the alternative H.  Of those two find the lowest and use it to create REC score.
import sys
import copy
import math
from scipy.stats.distributions import chi2

#**************************************************************************
def approxchi(mirInCanTypes):
    #summation of ranks
    approxchi = 0
    for el in mirInCanTypes:
            approxchi += math.log(el)
            #print( el, math.log(el) ,approxchi)

    print( "approxchi",approxchi *-2)
    return -2 * approxchi

#**************************************************************************
def getREC(pnul,pinvert):
    #Give the two pvals this give the rec score.
    #print( "getrec", pnul,pinvert

    if pinvert < pnul:
        print( "invRankPval smaller so Rec:")
        return math.log(2*pinvert,10)

    elif pinvert > pnul:
        print( "rankPval smaller so Rec:")
        return math.log(2*pnul,10) * -1

    elif pinvert == pnul:
        return 0
    else:
        return 0

#**************************************************************************

def main():
    mirdict = {}
    cancers = []
    Mir =[]
    findmir =[]
    with open(sys.argv[1]) as fin:
        cancers = fin.readline().split()
        findmir = fin.readline().split()
        for typ in cancers:
            mirdict[typ] = []

        Mir = [list(line.rstrip().split()) for line in fin]

    for key in mirdict:
        #makes a list if rvals etc. for each cancer
        mirdict[key] = [line for line in Mir if line[0] == key]
        # print( mirdict[key])

    for key in mirdict:
        n = len(mirdict[key])
        mirdict[key]= sorted(mirdict[key], key=lambda elem:elem[2] )
        # print( mirdict[key])

        for x,el in enumerate(mirdict[key]):
            el.append((float(x+1)/n) - (float(1)/(2*n)))
            el.append(float(n - (x+1) + 1)/n - (float(1)/(2*n)))

        mirdict[key] = [tuple(el) for el in mirdict[key]]
        # print( key,mirdict[key])
    MirNew = []

    for key in mirdict:
        MirNew += mirdict[key]
        #print( key,mirdict[key]



    for mir in findmir:
        apchi = 0
        #list of mir ranks
        mirInCanTypesTups = [tup for tup in MirNew if tup[1] == mir]

        mirInCanTypes = [tup[4] for tup in MirNew if tup[1] == mir]
        print( mir)
        #print( mirInCanTypesTups, '\n'
        df = len(mirInCanTypes)

        apchi = approxchi(mirInCanTypes)

        pval = chi2.sf(apchi, df*2)

        # print( "rankChiScore",apchi ,"df", df * 2, "rankpval",pval)

        #list of inverted mir ranks
        invertedMirInCanTypes = [tup[5] for tup in MirNew if tup[1] == mir]
        #print( invertedMirInCanTypes

        df = len(invertedMirInCanTypes)

        invapchi = approxchi(invertedMirInCanTypes)
        invpval = chi2.sf(invapchi,df*2)

        # print( "invRankChiScore",invapchi ,"df", df * 2, 'invertedrankpval', invpval)

        REC = getREC(pval,invpval)
        if REC < -3:
            print( "{0}\t{1}\tnumberOfCancers{2}\tGood Score\n".format(mir,REC,df))
        else:
            print( "{0}\t{1}\tnumberOfCancers{2}\n".format(mir,REC,df))

if __name__ == "__main__":
    main()







