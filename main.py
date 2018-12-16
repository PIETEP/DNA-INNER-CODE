#!/usr/bin/env python
import sys
import math
import numpy as np
from numpy.random import *
import factor_graph

def Readparam(filename):
    parameterlist=[]
    for line in open(filename,'r'):
        parameterlist.append(line)
    block_length=int(parameterlist[0])
    pa = float(parameterlist[1])
    pid = float(parameterlist[2])
    maxdrift = int(parameterlist[3])
    ps = float(parameterlist[4])
    nofcopy = int(parameterlist[5])
    Decordtimes =int(parameterlist[6])
    return block_length,pa,pid,maxdrift,ps,nofcopy,Decordtimes

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
    print(cordword)
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

def argmax_ndim(arg_array):
    return np.unravel_index(arg_array.argmax(), arg_array.shape)



def transmit(block_length,pa,pid,maxdrift,cordword,ps):
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
        if drift[i]==maxdrift:
            if rand() <= pid:
                drift.append(drift[i]-1)
            else:
                drift.append(drift[i])
        elif drift[i]==-maxdrift:
            if rand() <= pid:
                drift.append(drift[i]+1)
            else:
                drift.append(drift[i])
        else:
            if rand() <= pid:
                drift.append(drift[i]+1)
            elif rand() <= pid:
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


#sys.argv[1] is Channel Parameter file
block_length,pa,pid,maxdrift,ps,nofcopy,decordtimes=Readparam(sys.argv[1])
Allowable_error_count=int(nofcopy-1//2)
#sys,argv[2] is prior file
priorlist,period=Originalprior(sys.argv[2])
entropy(priorlist)
cordword=GetRndWord(block_length,period,priorlist)
#PriorCordword=SetPriorCordword(cordword,block_length,priorlist,period,nofcopy)

'''
downmessage=[[0 for i in range(block_length)] for j in range(nofcopy)]
leftmessage=[[0 for i in range(block_length)] for j in range(nofcopy)]
rightmessage=[[0 for i in range(block_length)] for j in range(nofcopy)]
'''

#print(upmessage)

CopyZ=[]
for copy_index in range(nofcopy):
    CopyZ.append(transmit(block_length, pa, pid, maxdrift, cordword, ps))
print(CopyZ)
NextPrior=np.tile(priorlist,(int(block_length/period),1))
#print('NextPrior',NextPrior)
for decode_index in range(decordtimes):
    for copy_index in range(nofcopy):
        #print(CopyZ[copy_index])

        normalized_out_probability=factor_graph.factor_graph(block_length,CopyZ[copy_index],maxdrift,NextPrior,period,pid,ps,priorlist)
        if copy_index!=0:
            out_probability_list=np.dstack([out_probability_list,normalized_out_probability])
        else:
            out_probability_list=normalized_out_probability
    #print(out_probability_list)

    #print(out_probability_list[0][0])
    #in case nofcopy=3
    for i in range(NextPrior.shape[0]):
        for j in range(NextPrior.shape[1]):
            NextPrior[i][j]=out_probability_list[i][j][0]*out_probability_list[i][j][1]*out_probability_list[i][j][2]+\
                            out_probability_list[i][j][0]*out_probability_list[i][j][1]*(1-out_probability_list[i][j][2])+\
                            out_probability_list[i][j][0]*(1-out_probability_list[i][j][1])*out_probability_list[i][j][2]+\
                            (1-out_probability_list[i][j][0])*out_probability_list[i][j][1]*out_probability_list[i][j][2]
    #print(NextPrior)
    sum_NextPrior=np.sum(NextPrior,axis=1)
    #print(sum_NextPrior)
    div=(np.tile(sum_NextPrior,(4,1))).T
    #print(div)
    NextPrior=NextPrior/div
    #print(NextPrior)

decode_word=[]
errorcount=0
for i in range(NextPrior.shape[0]):
    predicted_base=argmax_ndim(NextPrior[i])
    #print(predicted_base[0])
    decode_word.append(bin(predicted_base[0]))
    if decode_word[i]!=cordword[i]:
        errorcount+=1

print(decode_word)
print(errorcount/block_length)

