import numpy as np
cimport numpy as cnp

def factor_graph(int block_length,list Z,int maxdrift,int nofc,cnp.ndarray NextPrior,int period,double pd,double pi,double ps,list priorlist,cnp.ndarray ZProbTable):
    upmessage=NextPrior
    #print(upmessage,Z)
    #initial drift is zero
    rightmessage=SetUprightmessage(block_length,maxdrift,nofc,period,pd,pi,ps,upmessage,Z,ZProbTable)
    leftmessage=SetUpleftmessage(block_length,maxdrift,nofc,period,pd,pi,ps,upmessage,Z,ZProbTable)
    downmessage=SetUpDownmessage(block_length,maxdrift,period,pd,pi,ps,priorlist,Z,leftmessage,rightmessage,ZProbTable)
    out_proability=np.array(upmessage)*np.array(downmessage)
    sum_out_probability=np.sum(out_proability,axis=1)
    #print(out_proability)
    #print(sum_out_probability)
    div=(np.tile(sum_out_probability,(4,1))).T
    #print(div)
    normalize_out_probability=out_proability/div
    #print(normalize_out_probability)
    return normalize_out_probability

def SetUpleftmessage(int block_length,int maxdrift,int nofc,int period,double pd,double pi,double ps,cnp.ndarray upmessage,list Z,cnp.ndarray ZProbTable):
    cdef int i,j,k,di,nextdi,xi,firstdrift,seconddrift,thirddrift,firstnextdrift,secondnextdrift,thirdnextdrift,nofs
    cdef double ditonextdiprob,aboutZprobability,beforeleftmessage
    left=[[[[0 for thirddrift in range(2*maxdrift+1)] for seconddrift in range(2*maxdrift+1)] for firstdrift in range(2*maxdrift+1)] for j in range(block_length)]
    cdef cnp.ndarray leftmessage
    leftmessage=np.array(left,dtype=float)
    for firstdrift in range(2*maxdrift+1):
        for seconddrift in range(2*maxdrift+1):
            for thirddrift in range(2*maxdrift+1):
                leftmessage[block_length-1][firstdrift][seconddrift][thirddrift] = 1
    #print("deciding rightmessage")
    for i in range(block_length - 1,0,-1):
        if i < maxdrift:
            for firstdrift in range(-i-1, i + 2):
                for firstnextdrift in range(firstdrift - 1, firstdrift + 2):
                    if i + firstnextdrift > len(Z[0]) - 1 or firstnextdrift < -i or firstnextdrift >=i+1:
                        continue;
                    for seconddrift in range(-i-1, i + 2):
                        for secondnextdrift in range(seconddrift - 1, seconddrift + 2):
                            if i + secondnextdrift > len(Z[1]) - 1 or secondnextdrift < -i or secondnextdrift >=i+1:
                                continue;
                            for thirddrift in range(-i-1, i + 2):
                                for thirdnextdrift in range(thirddrift - 1, thirddrift + 2):
                                    if i + thirdnextdrift > len(Z[2]) - 1 or thirdnextdrift < -i or thirdnextdrift >=i+1:
                                        continue;
                                    beforeleftmessage= leftmessage[i][firstdrift + maxdrift][seconddrift +maxdrift][thirddrift+maxdrift]
                                    #print(thirdnextdrift+maxdrift)
                                    for xi in range(4):
                                        nofs=0
                                        for possibles in range(0,firstnextdrift-firstdrift+1):
                                            if bin(xi)!=Z[0][i+firstnextdrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,secondnextdrift-seconddrift+1):
                                            if bin(xi)!=Z[1][i+secondnextdrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,thirdnextdrift-thirddrift+1):
                                            if bin(xi)!=Z[2][i+thirdnextdrift-possibles]:
                                                nofs+=1
                                        #print(nofs)
                                        leftmessage[i-1][firstnextdrift + maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]+=upmessage[i][xi] *\
                                        ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1][seconddrift+maxdrift][secondnextdrift-seconddrift+1][thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*\
                                        beforeleftmessage
                                        #if i==0:
                                            #a=(ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1][seconddrift+maxdrift][secondnextdrift-seconddrift+1][thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*upmessage[i][xi])
                                            #print(a,firstdrift+maxdrift,firstnextdrift-firstdrift+1,seconddrift+maxdrift,secondnextdrift-seconddrift+1,thirddrift+maxdrift,thirdnextdrift-thirddrift+1,nofs) 
                                            #print(a*beforerightmessage,firstnextdrift+maxdrift,secondnextdrift+maxdrift,thirdnextdrift+maxdrift,rightmessage[1][firstnextdrift + maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]) 
        else:
            for firstdrift in range(-maxdrift, maxdrift+1):
                if i + firstdrift > len(Z[0]) - 1:
                    continue;
                for firstnextdrift in range(firstdrift - 1, firstdrift + 2):
                    if i + firstnextdrift > len(Z[0]) - 1 or firstnextdrift > maxdrift or firstnextdrift <-maxdrift:
                        continue;
                    for seconddrift in range(-maxdrift, maxdrift+1):
                        if i + seconddrift > len(Z[1]) - 1:
                            continue;
                        for secondnextdrift in range(seconddrift - 1, seconddrift + 2):
                            if i + secondnextdrift > len(Z[1]) - 1 or secondnextdrift > maxdrift or secondnextdrift <-maxdrift:
                                continue;
                            for thirddrift in range(-maxdrift, maxdrift+1):
                                if i + thirddrift > len(Z[2]) - 1:
                                    continue;
                                for thirdnextdrift in range(thirddrift - 1, thirddrift + 2):
                                    if i + thirdnextdrift > len(Z[2]) - 1 or thirdnextdrift > maxdrift or thirdnextdrift <-maxdrift:
                                        continue;
                                    beforeleftmessage= leftmessage[i][firstdrift + maxdrift][seconddrift +maxdrift][thirddrift+maxdrift]
                                    
                                    for xi in range(4):
                                        nofs=0
                                        for possibles in range(0,firstnextdrift-firstdrift+1):
                                            if bin(xi)!=Z[0][i+firstnextdrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,secondnextdrift-seconddrift+1):
                                            if bin(xi)!=Z[1][i+secondnextdrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,thirdnextdrift-thirddrift+1):
                                            if bin(xi)!=Z[2][i+thirdnextdrift-possibles]:
                                                nofs+=1

                                        leftmessage[i-1][firstnextdrift + maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]+=upmessage[i][xi] *\
                                        ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1]\
                                        [seconddrift+maxdrift][secondnextdrift-seconddrift+1]\
                                        [thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*\
                                        beforeleftmessage

        #normalization processi
        leftmessage[i-1]=sumone(leftmessage[i-1])
        #if i<5:
            #print(rightmessage[i][maxdrift][maxdrift])
    print(leftmessage)
    return leftmessage

def SetUprightmessage(int block_length,int maxdrift,int nofc,int period,double pd,double pi,double ps,cnp.ndarray upmessage,list Z, cnp.ndarray ZProbTable):
    cdef int i,j,k,di,nextdi,xi,firstdrift,seconddrift,thirddrift,firstnextdrift,secondnextdrift,thirdnextdrift,nofs
    cdef double ditonextdiprob,aboutZprobability
    right=[[[[0 for thirddrift in range(2*maxdrift+1)] for seconddrift in range(2*maxdrift+1)] for firstdrift in range(2*maxdrift+1)] for j in range(block_length)]
    cdef cnp.ndarray rightmessage
    rightmessage=np.array(right,dtype=float)
    rightmessage[0][maxdrift][maxdrift][maxdrift] = 1
    #print("deciding rightmessage")
    for i in range(0, block_length - 1):
        if i < maxdrift:
            for firstdrift in range(-i, i + 1):
                for firstnextdrift in range(firstdrift - 1, firstdrift + 2):
                    for seconddrift in range(-i, i + 1):
                        for secondnextdrift in range(seconddrift - 1, seconddrift + 2):
                            for thirddrift in range(-i, i + 1):
                                for thirdnextdrift in range(thirddrift - 1, thirddrift + 2):
                                    beforerightmessage= rightmessage[i][firstdrift + maxdrift][seconddrift +maxdrift][thirddrift+maxdrift]
                                    #print(thirdnextdrift+maxdrift)
                                    for xi in range(4):
                                        nofs=0
                                        for possibles in range(0,firstdrift-firstnextdrift+1):
                                            if bin(xi)!=Z[0][i+firstdrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,seconddrift-secondnextdrift+1):
                                            if bin(xi)!=Z[1][i+seconddrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,thirddrift-thirdnextdrift+1):
                                            if bin(xi)!=Z[2][i+thirddrift-possibles]:
                                                nofs+=1

                                        rightmessage[i + 1][firstnextdrift + maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]+=upmessage[i][xi] *\
                                        ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1][seconddrift+maxdrift][secondnextdrift-seconddrift+1][thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*\
                                        beforerightmessage
                                        #if i==0:
                                            #a=(ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1][seconddrift+maxdrift][secondnextdrift-seconddrift+1][thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*upmessage[i][xi])
                                            #print(a,firstdrift+maxdrift,firstnextdrift-firstdrift+1,seconddrift+maxdrift,secondnextdrift-seconddrift+1,thirddrift+maxdrift,thirdnextdrift-thirddrift+1,nofs) 
                                            #print(a*beforerightmessage,firstnextdrift+maxdrift,secondnextdrift+maxdrift,thirdnextdrift+maxdrift,rightmessage[1][firstnextdrift + maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]) 
        else:
            for firstdrift in range(-maxdrift, maxdrift+1):
                if i + firstdrift > len(Z[0]) - 1:
                    continue;
                for firstnextdrift in range(firstdrift - 1, firstdrift + 2):
                    if i + firstnextdrift > len(Z[0]) - 1 or firstnextdrift > maxdrift or firstnextdrift <-maxdrift:
                        continue;
                    for seconddrift in range(-maxdrift, maxdrift+1):
                        if i + seconddrift > len(Z[1]) - 1:
                            continue;
                        for secondnextdrift in range(seconddrift - 1, seconddrift + 2):
                            if i + secondnextdrift > len(Z[1]) - 1 or secondnextdrift > maxdrift or secondnextdrift <-maxdrift:
                                continue;
                            for thirddrift in range(-maxdrift, maxdrift+1):
                                if i + thirddrift > len(Z[2]) - 1:
                                    continue;
                                for thirdnextdrift in range(thirddrift - 1, thirddrift + 2):
                                    if i + thirdnextdrift > len(Z[2]) - 1 or thirdnextdrift > maxdrift or thirdnextdrift <-maxdrift:
                                        continue;
                                    beforerightmessage= rightmessage[i][firstdrift + maxdrift][seconddrift +maxdrift][thirddrift+maxdrift]
                                    #print(thirdnextdrift+maxdrift,rightmessage[i+1][firstnextdrift+maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift])
                                    for xi in range(4):
                                        nofs=0
                                        for possibles in range(0,firstdrift-firstnextdrift+1):
                                            if bin(xi)!=Z[0][i+firstdrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,seconddrift-secondnextdrift+1):
                                            if bin(xi)!=Z[1][i+seconddrift-possibles]:
                                                nofs+=1
                                        for possibles in range(0,thirddrift-thirdnextdrift+1):
                                            if bin(xi)!=Z[2][i+thirddrift-possibles]:
                                                nofs+=1

                                        rightmessage[i + 1][firstnextdrift + maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]+=upmessage[i][xi] *\
                                        ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1]\
                                        [seconddrift+maxdrift][secondnextdrift-seconddrift+1]\
                                        [thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*\
                                        beforerightmessage

        #normalization processi
        rightmessage[i+1]=sumone(rightmessage[i+1])
        if i%10==1:
            print(rightmessage[i][maxdrift][maxdrift])
    print(rightmessage)
    return rightmessage

def SetUpDownmessage(int block_length,int maxdrift,int period,double pd,double pi,double ps,list priorlist,list Z,cnp.ndarray leftmessage,cnp.ndarray rightmessage,cnp.ndarray ZProbTable):
    cdef int i,j,k,limit,di,beforedi,xi,firstdrift,firstnextdrift,seconddrift,secondnextdrift,thirddrift,thirdnextdrift
    cdef list down=[[0 for i in range(4)] for j in range(block_length)]
    cdef cnp.ndarray downmessage
    downmessage=np.array(down,dtype=float)
    for i in range(block_length):
        if i<maxdrift:
            limit=i
        else:
            limit=maxdrift
        for firstdrift in range(-limit,limit+1):
            if i+firstdrift>len(Z[0])-1:
                continue;
            for firstnextdrift in range(firstdrift-1,firstdrift+2):
                if i+firstnextdrift>len(Z[0])-1 or firstnextdrift>maxdrift or firstnextdrift<-maxdrift:
                    continue;
                for seconddrift in range(-limit,limit+1):
                    if i+seconddrift>len(Z[1])-1:
                        continue;
                    for secondnextdrift in range(seconddrift-1,seconddrift+2):
                        if i+secondnextdrift>len(Z[1])-1 or secondnextdrift>maxdrift or secondnextdrift<-maxdrift:
                            continue;
                        for thirddrift in range(-limit,limit+1):
                            if i+thirddrift>len(Z[2])-1:
                                continue;
                            for thirdnextdrift in range(thirddrift-1,thirddrift+2):
                                if i+thirdnextdrift>len(Z[2])-1 or thirdnextdrift>maxdrift or thirdnextdrift<-maxdrift:
                                    continue;

                                for xi in range(4):
                                    nofs=0
                                    for possibles in range(0,firstdrift-firstnextdrift+1):
                                        if bin(xi)!=Z[0][i+firstdrift-possibles]:
                                            nofs+=1
                                    for possibles in range(0,seconddrift-secondnextdrift+1):
                                        if bin(xi)!=Z[1][i+seconddrift-possibles]:
                                            nofs+=1
                                    for possibles in range(0,thirddrift-thirdnextdrift+1):
                                        if bin(xi)!=Z[2][i+thirddrift-possibles]:
                                            nofs+=1
                                    #print(nofs,downmessage[i][xi])
                                    downmessage[i][xi]+=ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1]\
                                                        [seconddrift+maxdrift][secondnextdrift-seconddrift+1]\
                                                        [thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*\
                                                        rightmessage[i][firstdrift+maxdrift][seconddrift+maxdrift][thirddrift+maxdrift]*\
                                                        leftmessage[i][firstnextdrift+maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]
                                    #print(i, ditonextdiprob,aboutZprobability,rightmessage[i][di+maxdrift],leftmessage[i][nextdi+maxdrift],downmessage[i][xi])
        downmessage[i] = sumone(downmessage[i])
    print(downmessage)
    return downmessage

def sumone(x,axis=None):
    a=np.sum(x)
    result=x/a
    return result
