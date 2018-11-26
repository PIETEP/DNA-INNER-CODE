#!/usr/bin/env python
import sys
import math
import numpy as np
from numpy.random import *

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
    print("SetPriorCordword",PriorCordword)

    return PriorCordword


'''
def Setupmessage():
    if decode_index>=1:
        x=5
    else:
        return PriorCordword
'''

def SetUpleftmessage(block_length,maxdrift,period,pid,ps,priorlist,Z):
    leftmessage = [[0 for i in range(2 * maxdrift + 1)] for j in range(block_length)]
    for k in range(2*maxdrift+1):
        leftmessage[block_length-1][k]=1
    print("deciding leftmessage")
    for i in range(block_length - 1,0,-1):
        for xi in range(4):
            if i < 3:
                for di in range(-1-i, i + 2):
                    if i + di > len(Z) - 1:
                        continue;
                    for beforedi in range(di-1, di + 2):
                        if i + beforedi > len(Z) - 1 or beforedi < -i or beforedi >=i+1:
                            continue;
                        if beforedi == di:
                            ditonextdiprob = 1 - pid
                            if bin(xi) == Z[i+beforedi]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                        elif beforedi == di - 1:
                            ditonextdiprob = pid
                            aboutZprobability = 1
                        else:
                            ditonextdiprob = pid
                            if bin(xi) == Z[i + beforedi]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                            if bin(xi) == Z[i + beforedi - 1]:
                                aboutZprobability *= (1 - ps)
                            else:
                                aboutZprobability *= (ps / 3)
                        #print(i,beforedi)
                        leftmessage[i - 1][beforedi+maxdrift] += priorlist[i % period][xi] * ditonextdiprob * aboutZprobability * \
                                               leftmessage[i][beforedi+maxdrift]
            else:
                for di in range(-maxdrift, maxdrift+1):
                    if i + di > len(Z) - 1:
                        continue;
                    for beforedi in range(di-1, di + 2):
                        if i + beforedi > len(Z) - 1 or beforedi < -maxdrift or beforedi > maxdrift:
                            continue;
                        if beforedi == di:
                            ditonextdiprob = 1 - pid
                            if bin(xi) == Z[i+beforedi]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                        elif beforedi == di - 1:
                            ditonextdiprob = pid
                            aboutZprobability = 1
                        else:
                            ditonextdiprob = pid
                            if bin(xi) == Z[i + beforedi]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                            if bin(xi) == Z[i + beforedi - 1]:
                                aboutZprobability *= (1 - ps)
                            else:
                                aboutZprobability *= (ps / 3)
                        #print(beforedi)
                        leftmessage[i - 1][beforedi+maxdrift] += priorlist[i % period][xi] * ditonextdiprob * aboutZprobability * \
                                               leftmessage[i][di+maxdrift]
                        #print(i,di,beforedi,leftmessage[i-1][beforedi+maxdrift],leftmessage[i][beforedi+maxdrift])
    #print(leftmessage)
    return leftmessage


def SetUprightmessage(block_length,maxdrift,period,pid,ps,priorlist,Z):
    rightmessage=[[0 for i in range(2*maxdrift+1)] for j in range(block_length)]
    rightmessage[0][maxdrift] = 1
    print("deciding rightmessage")
    for i in range(0, block_length - 1):
        for xi in range(4):
            if i < maxdrift:
                for di in range(-i, i + 1):
                    if i + di > len(Z) - 1:
                        continue;
                    for nextdi in range(di - 1, di + 2):
                        if i + nextdi > len(Z) - 1:
                            continue;
                        if nextdi == di:
                            ditonextdiprob = 1 - pid
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3

                        elif nextdi == di - 1:
                            ditonextdiprob = pid
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                            if bin(xi) == Z[i + di - 1]:
                                aboutZprobability *= (1 - ps)
                            else:
                                aboutZprobability *= (ps / 3)
                        else:
                            ditonextdiprob = pid
                            aboutZprobability = 1
                        rightmessage[i + 1][nextdi + maxdrift] += priorlist[i % period][xi] * ditonextdiprob * aboutZprobability * \
                                                                 rightmessage[i][di + maxdrift]
            else:
                for di in range(-maxdrift, maxdrift+1):
                    if i + di > len(Z) - 1:
                        continue;
                    for nextdi in range(di - 1, di + 2):
                        if i + nextdi > len(Z) - 1 or nextdi > maxdrift or nextdi <-maxdrift:
                            continue;
                        if nextdi == di:
                            ditonextdiprob = 1 - pid
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                        elif nextdi == di - 1:
                            ditonextdiprob = pid
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                            if bin(xi) == Z[i + di - 1]:
                                aboutZprobability *= (1 - ps)
                            else:
                                aboutZprobability *= (ps / 3)
                        else:
                            ditonextdiprob = pid
                            aboutZprobability = 1
                        #print(rightmessage[i + 1], ditonextdiprob, aboutZprobability)
                        rightmessage[i + 1][nextdi + maxdrift] += priorlist[i % period][xi] * ditonextdiprob * aboutZprobability * \
                                                                 rightmessage[i][di + maxdrift]

    #print(rightmessage)
    return rightmessage

def SetUpDownmessage(block_length,maxdrift,period,pid,ps,priorlist,Z,leftmessage,rightmessage):
    downmessage=[[0 for i in range(4)] for j in range(block_length)]
    for i in range(block_length):
        if i<maxdrift:
            for di in range(-i,i+1):
                if i+di>len(Z)-1:
                    continue;
                for nextdi in range(di-1,di+2):
                    if i+nextdi>len(Z)-1:
                        continue;
                    elif di==nextdi:
                        ditonextdiprob=1-pid
                    else:
                        ditonextdiprob=pid
                    for xi in range(4):
                        if di==nextdi:
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                        elif nextdi==di-1:
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                            if bin(xi) == Z[i + di - 1]:
                                aboutZprobability *= (1 - ps)
                            else:
                                aboutZprobability *= (ps / 3)
                        else:
                            aboutZprobability = 1

                        downmessage[i][xi]+=ditonextdiprob*aboutZprobability*rightmessage[i][di+maxdrift]*leftmessage[i][nextdi+maxdrift]
                        #print(i, ditonextdiprob,aboutZprobability,rightmessage[i][di+maxdrift],leftmessage[i][nextdi+maxdrift],downmessage[i][xi])
        else:
            for di in range(-maxdrift,maxdrift+1):
                if i+di>len(Z)-1:
                    continue;
                for nextdi in range(di-1,di+2):
                    if i+nextdi>len(Z)-1 or nextdi>maxdrift or nextdi<-maxdrift:
                        continue;
                    elif di==nextdi:
                        ditonextdiprob=1-pid
                    else:
                        ditonextdiprob=pid
                    for xi in range(4):
                        if di==nextdi:
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                        elif nextdi==di-1:
                            if bin(xi) == Z[i + di]:
                                aboutZprobability = 1 - ps
                            else:
                                aboutZprobability = ps / 3
                            if bin(xi) == Z[i + di - 1]:
                                aboutZprobability *= (1 - ps)
                            else:
                                aboutZprobability *= (ps / 3)
                        else:
                            aboutZprobability = 1
                        downmessage[i][xi]+=ditonextdiprob*aboutZprobability*rightmessage[i][di+maxdrift]*leftmessage[i][nextdi+maxdrift]
    print(downmessage)
    return downmessage






def factor_graph(block_length,Z,copy_index,decode_index,maxdrift,period,pid,ps,priorlist):
    if decode_index == 0:
        upmessage = np.tile(priorlist,(int(block_length/period),1))
    else:
        upmessage = 5#NextPrior[copy_index]
    print(upmessage,Z)
    #initial drift is zero
    rightmessage=SetUprightmessage(block_length,maxdrift,period,pid,ps,priorlist,Z)
    leftmessage=SetUpleftmessage(block_length,maxdrift,period,pid,ps,priorlist,Z)
    downmessage=SetUpDownmessage(block_length,maxdrift,period,pid,ps,priorlist,Z,leftmessage,rightmessage)
    #NextPrior=



#sys.argv[1] is Channel Parameter file
block_length,pa,pid,maxdrift,ps,nofcopy,decordtimes=Readparam(sys.argv[1])
#sys,argv[2] is prior file
priorlist,period=Originalprior(sys.argv[2])
entropy(priorlist)
cordword=GetRndWord(block_length,period,priorlist)
PriorCordword=SetPriorCordword(cordword,block_length,priorlist,period,nofcopy)

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
for decode_index in range(decordtimes):
    for copy_index in range(nofcopy):
        #print(CopyZ[copy_index])


        factor_graph(block_length,CopyZ[copy_index],copy_index,decode_index,maxdrift,period,pid,ps,priorlist)

