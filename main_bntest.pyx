#!/usr/bin/env python
import sys
import math
import numpy as np
cimport numpy as cnp
from numpy.random import *
import decode
import factor_graph
import joint_graph
import setdriftprobtable

def Readparam(filename):
    parameterlist=[]
    for line in open(filename,'r'):
        parameterlist.append(line)
    block_length=int(parameterlist[0])
    pa = float(parameterlist[1])
    pi = float(parameterlist[2])
    pd = float(parameterlist[2])
    maxdrift = int(parameterlist[3])
    ps = float(parameterlist[4])
    nofcopy = int(parameterlist[5])
    Decordtimes =int(parameterlist[6])
    return block_length,pa,pi,pd,maxdrift,ps,nofcopy,Decordtimes

    #print(parameterlist)

def Originalprior(filename):
    priorlist=[]
    for line in open(filename,'r'):
        line = line.rstrip('\r\n')
        items = line.split(',')
        numbers = [float(item) for item in items] #str->float
        priorlist.append(numbers)
    period=len(priorlist)
    #print(period)
    return priorlist,period

def entropy(priorlist):
    entropy=[]
    for H in priorlist:
        temp=0
        for i in range(4):
            if H[i]==0.0:
                continue
            else:
                temp-=(H[i]*math.log2(H[i]))
        entropy.append(temp)
    #print(entropy)


def GetRndWord(block_length,period,priorlist):
    cordword=[]
    for i in range(block_length):
        flag=0
        count=0
        pool=0
        temp=rand()
        #Determine base by generated random number
        while(flag==0):
            #print(i%period,count)
            pool+=priorlist[i % period][count]
            if temp <= pool:
                cordword.append(bin(count))
                flag=1
            count+=1
    #print(cordword)
    return cordword

def substitution(enki,probability):
    temp = rand()
    Nucleotide = [bin(0), bin(1), bin(2), bin(3)]
    Nucleotide.remove(enki)
    if temp <= (1 / probability):
        return Nucleotide[0]
    elif temp <= (2 / probability):
        return Nucleotide[1]
    else:
        return Nucleotide[2]



def transmit(block_length,pa,pd,pi,maxdrift,cordword,ps):
    drift=[int(0)]
    state=[]
    Zprime=[]
    Z=[]
    for i in range(block_length+maxdrift):
        #state
        if rand()<=pa:
            state.append(1)
        else:
            state.append(0)
        #drift
        temp=rand()
        if drift[i]==maxdrift:
            if temp <= pi:
                drift.append(drift[i]-1)
            else:
                drift.append(drift[i])
        elif drift[i]==-maxdrift:
            if temp <= pd:
                drift.append(drift[i]+1)
            else:
                drift.append(drift[i])
        else:
            if temp <= pd:
                drift.append(drift[i]+1)
            elif temp >= 1-pi:
                drift.append(drift[i]-1)
            else:
                drift.append(drift[i])
        #Zprime
        if i+drift[i]>=block_length:
            break;
        else:
            Zprime.append(cordword[i+drift[i]])
            #Z
            if state[i]==1:
                Z.append(substitution(Zprime[i],3))
            else:
                if rand()>ps:
                    Z.append(Zprime[i])
                else:
                    Z.append(substitution(Zprime[i],3))


    #print(drift,state,Zprime,Z)
    return Z

def SetPriorCordword(cordword,block_length,priorlist,period,nofcopy):
    PriorCordword=[]
    for i in range(block_length):
        for j in range(4):
            if cordword[i]==bin(j):
                PriorCordword.append(priorlist[i%period][j])
    #print("SetPriorCordword",PriorCordword)

    return PriorCordword


'''
def Setupmessage():
    if decode_index>=1:
        x=5
    else:
        return PriorCordword
'''


def simulate(a,b):
    #sys.argv[1] is Channel Parameter file
    block_length,pa,pi,pd,maxdrift,ps,nofcopy,decordtimes=Readparam(a)
    Allowable_error_count=int(nofcopy-1//2)
    #sys,argv[2] is prior file
    priorlist,period=Originalprior(b)
    entropy(priorlist)
    #PriorCordword=SetPriorCordword(cordword,block_length,priorlist,period,nofcopy)


    #print(upmessage)

    wordcount=0
    worderrorcount=0
    rawerrorcount=0
    decodeerrorcount=0
    decode_word_error_count=0

    #print(setdriftprobtable.setprobtable(maxdrift,pd,pi,ps))
    ZProbTable=setdriftprobtable.setprobtable(maxdrift,pd,pi,ps)
    ZProbTable=np.array(ZProbTable,dtype=float)
    while(wordcount<1000):
        wordcount+=1
        cordword=GetRndWord(block_length,period,priorlist)
        CopyZ = []
        for copy_index in range(nofcopy):
            CopyZ.append(transmit(block_length, pa, pd,pi,maxdrift, cordword, ps))
            flag=0
            for i in range(min(len(cordword), len(CopyZ[copy_index]))):
                if cordword[i] != CopyZ[copy_index][i]:
                    rawerrorcount += 1
                    if flag==0:
                        flag=1
                        worderrorcount+=1
        print(CopyZ)
        NextPrior=np.tile(priorlist,(int(block_length/period),1))
        #print('initialPrior',NextPrior,NextPrior.shape)
                #print(decode_word)
        
        #decoderrorcount_plus,decode_word_error_count_plus=decode(block_length,CopyZ,cordword,decordtimes,maxdrift,nofcopy,NextPrior,period,pid,ps,priorlist)
        decodeerrorcount_plus,decode_word_error_count_plus=decode.joint_decode(block_length,cordword,CopyZ,maxdrift,nofcopy,NextPrior,period,pd,pi,ps,priorlist,ZProbTable)
        
        decodeerrorcount+=decodeerrorcount_plus
        decode_word_error_count+=decode_word_error_count_plus
        print(wordcount)
        print(worderrorcount)
        print(decode_word_error_count)

    print('Word Count:%d' % (wordcount))
    print('Raw Word Error count:%d' % (worderrorcount))
    print('decode Word error count:%d' % (decode_word_error_count))
    print('Raw Word Error Rate:%f' % (worderrorcount/(wordcount*nofcopy)))
    print('raw Symbol Error Rate:%f' % (rawerrorcount/(block_length*nofcopy*wordcount)))
    print('decode Word Error Rate:%f' %(decode_word_error_count/wordcount))
    print('decode Symbol Error Rate:%f' % (decodeerrorcount/(block_length*wordcount)))

