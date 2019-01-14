import numpy as np
cimport  numpy as cnp

def factor_graph(int block_length,list Z,int maxdrift,cnp.ndarray NextPrior,int period,double pid,double ps,list priorlist):
    upmessage=NextPrior
    #print(upmessage,Z)
    #initial drift is zero
    rightmessage=SetUprightmessage(block_length,maxdrift,period,pid,ps,upmessage,Z)
    leftmessage=SetUpleftmessage(block_length,maxdrift,period,pid,ps,upmessage,Z)
    downmessage=SetUpDownmessage(block_length,maxdrift,period,pid,ps,priorlist,Z,leftmessage,rightmessage)
    out_proability=np.array(upmessage)*np.array(downmessage)
    sum_out_probability=np.sum(out_proability,axis=1)
    #print(out_proability)
    #print(sum_out_probability)
    div=(np.tile(sum_out_probability,(4,1))).T
    #print(div)
    normalize_out_probability=out_proability/div
    #print(normalize_out_probability)
    return normalize_out_probability

def SetUpleftmessage(int block_length,int maxdrift,int period,double pid,double ps,cnp.ndarray upmessage,list Z):
    cdef int i,j,k,di,beforedi,xi
    cdef double ditonextdiprob,aboutZprobability
    cdef list leftmessage = [[0 for i in range(2 * maxdrift + 1)] for j in range(block_length)]


    for k in range(2*maxdrift+1):
        leftmessage[block_length-1][k]=1
    #print("deciding leftmessage")
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
                        leftmessage[i - 1][beforedi+maxdrift] += upmessage[i][xi] * ditonextdiprob * aboutZprobability * \
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
                        leftmessage[i - 1][beforedi+maxdrift] += upmessage[i][xi] * ditonextdiprob * aboutZprobability * \
                                               leftmessage[i][di+maxdrift]
                        #print(i,di,beforedi,leftmessage[i-1][beforedi+maxdrift],leftmessage[i][beforedi+maxdrift])
        # normalization process
        leftmessage_array = np.array(leftmessage[i - 1])
        leftmessage_array_normalize = sumone(leftmessage_array)
        leftmessage[i - 1] = leftmessage_array_normalize.tolist()
    #print(leftmessage)
    return leftmessage


def SetUprightmessage(int block_length,int maxdrift,int period,double pid,double ps,cnp.ndarray upmessage,list Z):
    cdef int i,j,k,di,nextdi,xi
    cdef double ditonextdiprob,aboutZprobability
    cdef list rightmessage=[[0 for i in range(2*maxdrift+1)] for j in range(block_length)]
    rightmessage[0][maxdrift] = 1
    #print("deciding rightmessage")
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
                        rightmessage[i + 1][nextdi + maxdrift] += upmessage[i][xi] * ditonextdiprob * aboutZprobability * \
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
                        rightmessage[i + 1][nextdi + maxdrift] += upmessage[i][xi] * ditonextdiprob * aboutZprobability * \
                                                                 rightmessage[i][di + maxdrift]
        #normalization process
        rightmessage_array=np.array(rightmessage[i+1])
        rightmessage_array_normalize=sumone(rightmessage_array)
        rightmessage[i+1]=rightmessage_array_normalize.tolist()
    #print(rightmessage)
    return rightmessage

def SetUpDownmessage(int block_length,int maxdrift,int period,double pid,double ps,list priorlist,list Z,list leftmessage,list rightmessage):
    cdef int i,j,k,di,beforedi,xi
    cdef double ditonextdiprob,aboutZprobability
    cdef list downmessage=[[0 for i in range(4)] for j in range(block_length)]
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
        downmessage_array = np.array(downmessage[i])
        downmessage_array_normalize = sumone(downmessage_array)
        downmessage[i] = downmessage_array_normalize.tolist()
    #print(downmessage)
    return downmessage

def sumone(x,axis=None):
    a=np.sum(x)
    result=x/a
    return result