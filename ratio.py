#python ratio.py L1HS.cancer.txt L1HS.normal.txt miRnaCountsPerMillion.txt
normal.ReadCountPerMillion patientNCancer.csv> output

import sys
import math

from copy import deepcopy
from scipy import stats
from itertools import chain
from Bycancerfun import sortPatientbyCancer

import matplotlib
import matplotlib.pyplot as plt

matplotlib.style.use('ggplot')

#**************************************************************************
class miRna:
    def __init__(self,can,name,rval,pval,slope,intercept,Xaxis,Yaxis):
        self.name = name
        self.canType = can
        self.rval = rval
        self.pval = pval
        self.slope = slope
        self.intercept = intercept
        self.Xaxis = Xaxis
        self.Yaxis = Yaxis

#**************************************************************************
def handleZeroTup(line1,line2):
    '''expects 2 lists of tuples.  since it's the ratio the (a:b = (a,b)).
    Since we are make 2 axis there will be 2 lists.  If there is a zero in
    thefirst tup the whole line x,y point will be invalidated'''

    lineXY= list(zip(line1,line2))

    lineXY= [tups for tups in lineXY if not (0 in list(chain(*tups))) ]

    return lineXY
#**************************************************************************
def makeMirDict(mirTup,keys,TypesDictionary):
    curMirDictionary = dict((k,[]) for k in keys) #make dictionary

    for can in keys:
        curMirDictionary[can] = [mirTup[ind] for ind in TypesDictionary[can]] #make each dictionary
        #by getting the index lists you created.

    return curMirDictionary

#**************************************************************************
def plotReg(name,can,rval,pval,Xaxis,Yaxis,slope,intercept,col):
        title = "{0},pVal:{2:10.4f}".format(name)
        plt.title(title)# add title to graph
        plt.xlabel(name + "  reads per million MiRNA mapped")# add X label
        plt.ylabel("Transposon Expression") # Y label

        x1 = min(Xaxis)
        x2 = max(Xaxis)
        y1 = x1 * slope + intercept
        y2 = x2 * slope + intercept

        plt.scatter(Xaxis,Yaxis, c= col,edgecolor = "black")
        plt.plot([x1, x2], [y1, y2], c = col) # plot transposons and MiRna counts
        plt.savefig('graphs3/'+ title + '.png', format='png', dpi=500) # save to file
        plt.close()
        return None

#**************************************************************************
#a sort by with key function will be useful
def main():
    canTrans = sys.argv[1]
    normTrans = sys.argv[2]
    canMir = sys.argv[3]
    normMir = sys.argv[4]
    patientCan = sys.argv[5]

    patients = []
    cancertran = []
    normaltran = []
    indexList = []
    listMiRna =[]

    TypesDictionary = sortPatientbyCancer(patientCan)
    #This funtion is the most important part of the program
    #For every cancer it goes through the list, marks indexes with that cancertype
    #When buiding lists in the future, it looks up which patients with this dictionary
    keys = TypesDictionary.keys()
    #print("\t".join(keys))


    with open(normTrans) as norms, open(canTrans) as cans:
        patients = cans.readline().split()[1:]
        norms.readline() #throw away

        cancertran = cans.readline().split()[1:]
        normaltran = norms.readline().split()[1:]


        cancertran = list(map(float, cancertran))
        normaltran = list(map(float, normaltran))

    with open(normMir) as normsM, open(canMir) as cansM:
        mirPatients = cansM.readline().split()[1:]
        normsM.readline() #throw away


        for patient in mirPatients: #this is to sort by transPatinet(patients) by their
            #position in mirpatient
            ind = patients.index(patient)#find where trans patient corresponds to mir patient
            indexList.append(ind)

        cancertran = [cancertran[ind] for ind in indexList]#puts element in correct place
        normaltran = [normaltran[ind] for ind in indexList]#puts element in correct place

        tranTup = list(zip(cancertran,normaltran))

        TranByCanDictionary = makeMirDict(tranTup,keys,TypesDictionary)
        #print( TranByCanDictionary['THYM']

        for line in cansM:
            normLine = normsM.readline().split()[1:]
            normLine = list(map(float,normLine))

            canLine = line.split()
            name = canLine[0]
            canLine = list(map(float,canLine[1:]))

            mirTup = list(zip(canLine,normLine))
            curMirDictionary = makeMirDict(mirTup,keys,TypesDictionary)
            #Creates dictionary for cancer type build each with list of indexes
            noZeroDictionary = dict((k,[]) for k in keys)

            tempTranByCanDictionary = deepcopy(TranByCanDictionary)

            for can in TypesDictionary:
                noZeroDictionary[can] =  handleZeroTup(curMirDictionary[can],tempTranByCanDictionary[can])
                #So the function above zips the lists two, chains the resulting tuples,  if ther is a zero it's deleted

            for can in TypesDictionary:
                if len(noZeroDictionary[can]) > 0:
                    curMirDictionary[can], tempTranByCanDictionary[can] = list(zip(*noZeroDictionary[can]))#unzip
                    curMirDictionary[can] = list(curMirDictionary[can]) # get list of tuples of mir's (cancer,normal)
                    tempTranByCanDictionary[can] = list(tempTranByCanDictionary[can]) #get list of tuples of trans(cancernormal
                    curMirDictionary[can] = [math.log(cancer/normal,2) for cancer,normal in list(curMirDictionary[can])] #log2(ratio)
                    tempTranByCanDictionary[can] = [math.log(cancer/normal,2) for cancer,normal in list(tempTranByCanDictionary[can])]
                    slope,intercept,rval,pval,_ = stats.linregress(curMirDictionary[can],tempTranByCanDictionary[can])


                    print( "{0}\t{1}\t{2:0.9f}\t{4:0.9f}\tslope:{5:0.9f}\tintercept{6:0.9f}\t{3}".format(can,name,pval,len(curMirDictionary[can]), rval,slope,intercept))

if __name__ == "__main__":
    main()
