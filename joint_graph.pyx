import numpy as np
cimport numpy as cnp

def factor_graph(int block_length,list drift,list Z,int maxdrift,int nofc,cnp.ndarray NextPrior,int period,double pd,double pi,double ps,cnp.ndarray ZProbTable):
    upmessage=NextPrior
    #print(upmessage,Z)
    #initial drift is zero
    rightmessage=SetUprightmessage(block_length,drift,maxdrift,nofc,period,pd,pi,ps,upmessage,Z,ZProbTable)
    leftmessage=SetUpleftmessage(block_length,drift,maxdrift,nofc,period,pd,pi,ps,upmessage,Z,ZProbTable)
    downmessage=SetUpDownmessage(block_length,maxdrift,period,pd,pi,ps,Z,leftmessage,rightmessage,ZProbTable)
    out_proability=np.array(upmessage)*np.array(downmessage)
    sum_out_probability=np.sum(out_proability,axis=1)
    #print(out_proability)
    #print(sum_out_probability)
    div=(np.tile(sum_out_probability,(4,1))).T
    #print(div)
    normalize_out_probability=out_proability/div
    #print(normalize_out_probability)
    return normalize_out_probability

def SetUpleftmessage(int block_length,list drift,int maxdrift,int nofc,int period,double pd,double pi,double ps,cnp.ndarray upmessage,list Z,cnp.ndarray ZProbTable):
    cdef int i,xi,firstdrift,seconddrift,thirddrift,firstbeforedrift,secondbeforedrift,thirdbeforedrift,nofs,Zindex
    cdef double beforeleftmessage
    cdef cnp.ndarray leftmessage
    messagesize=2*maxdrift+1
    leftmessage=np.zeros((block_length,messagesize,messagesize,messagesize),dtype=float)
    index=[]
    for i,d in enumerate(Z):
        index.append(block_length-len(Z[i])+maxdrift)
    print(index)
    indexmin=min(index)
    leftmessage[block_length-1][index[0]][index[1]][index[2]] = 1
    print("deciding leftmessage")
    for i in range(block_length - 1,0,-1):
        if i < maxdrift:
            for firstdrift in range(-i-1, i + 2):
                for firstbeforedrift in range(firstdrift - 1, firstdrift + 2):
                    if i - firstbeforedrift > len(Z[0]) - 1 or firstbeforedrift < -i or firstbeforedrift >i:
                        continue;
                    for seconddrift in range(-i-1, i + 2):
                        for secondbeforedrift in range(seconddrift - 1, seconddrift + 2):
                            if i - secondbeforedrift > len(Z[1]) - 1 or secondbeforedrift < -i or secondbeforedrift >i:
                                continue;
                            for thirddrift in range(-i-1, i + 2):
                                for thirdbeforedrift in range(thirddrift - 1, thirddrift + 2):
                                    if i - thirdbeforedrift > len(Z[2]) - 1 or thirdbeforedrift < -i or thirdbeforedrift >i:
                                        continue;
                                    beforeleftmessage= leftmessage[i][firstdrift + maxdrift][seconddrift +maxdrift][thirddrift+maxdrift]
                                    #print(thirdbeforedrift+maxdrift)
                                    nowdrift=np.array([firstdrift,seconddrift,thirddrift])
                                    beforedrift=np.array([firstbeforedrift,secondbeforedrift,thirdbeforedrift])
                                    if indexmin==0:
                                        print(beforeleftmessage,nowdrift,beforedrift)
                                    for xi in range(4):
                                        nofs=0
                                        for Zindex in range(3):
                                            for possibles in range(beforedrift[Zindex]-nowdrift[Zindex]+1):
                                                if len(Z[Zindex])-1<i-beforedrift[Zindex]+possibles:
                                                    continue;
                                                if str(xi)!=Z[Zindex][i-beforedrift[Zindex]+possibles]:
                                                    nofs+=1
                                        
                                        #print(nofs)
                                        leftmessage[i-1][firstbeforedrift + maxdrift][secondbeforedrift+maxdrift][thirdbeforedrift+maxdrift]+=upmessage[i][xi] *\
                                        ZProbTable[firstbeforedrift+maxdrift][firstdrift-firstbeforedrift+1][secondbeforedrift+maxdrift][seconddrift-secondbeforedrift+1][thirdbeforedrift+maxdrift][thirddrift-thirdbeforedrift+1][nofs]*\
                                        beforeleftmessage
                                        #if i==0:
                                            #a=(ZProbTable[firstdrift+maxdrift][firstbeforedrift-firstdrift+1][seconddrift+maxdrift][secondbeforedrift-seconddrift+1][thirddrift+maxdrift][thirdbeforedrift-thirddrift+1][nofs]*upmessage[i][xi])
                                            #print(a,firstdrift+maxdrift,firstbeforedrift-firstdrift+1,seconddrift+maxdrift,secondbeforedrift-seconddrift+1,thirddrift+maxdrift,thirdbeforedrift-thirddrift+1,nofs) 
                                            #print(a*beforerightmessage,firstbeforedrift+maxdrift,secondbeforedrift+maxdrift,thirdbeforedrift+maxdrift,rightmessage[1][firstbeforedrift + maxdrift][secondbeforedrift+maxdrift][thirdbeforedrift+maxdrift]) 
        else:
            for firstdrift in range(-maxdrift, maxdrift+1):
                if i - firstdrift > len(Z[0]) - 1:
                    continue;
                for firstbeforedrift in range(firstdrift - 1, firstdrift + 2):
                    if i - firstbeforedrift > len(Z[0]) - 1 or firstbeforedrift > maxdrift or firstbeforedrift <-maxdrift:
                        continue;
                    for seconddrift in range(-maxdrift, maxdrift+1):
                        if i - seconddrift > len(Z[1]) - 1:
                            continue;
                        for secondbeforedrift in range(seconddrift - 1, seconddrift + 2):
                            if i - secondbeforedrift > len(Z[1]) - 1 or secondbeforedrift > maxdrift or secondbeforedrift <-maxdrift:
                                continue;
                            for thirddrift in range(-maxdrift, maxdrift+1):
                                if i - thirddrift > len(Z[2]) - 1:
                                    continue;
                                for thirdbeforedrift in range(thirddrift - 1, thirddrift + 2):
                                    if i - thirdbeforedrift > len(Z[2]) - 1 or thirdbeforedrift > maxdrift or thirdbeforedrift <-maxdrift:
                                        continue;
                                    beforeleftmessage= leftmessage[i][firstdrift + maxdrift][seconddrift +maxdrift][thirddrift+maxdrift]
                                    nowdrift=np.array([firstdrift,seconddrift,thirddrift])
                                    beforedrift=np.array([firstbeforedrift,secondbeforedrift,thirdbeforedrift])
                                    #if indexmin==0 and beforeleftmessage!=0:
                                        #print(beforeleftmessage,nowdrift,beforedrift)
                                    for xi in range(4):
                                        nofs=0
                                        for Zindex in range(3):
                                            for possibles in range(beforedrift[Zindex]-nowdrift[Zindex]+1):
                                                if len(Z[Zindex])-1<i-beforedrift[Zindex]+possibles:
                                                    continue;
                                                if str(xi)!=Z[Zindex][i-beforedrift[Zindex]+possibles]:
                                                    nofs+=1
                                        leftmessage[i-1][firstbeforedrift + maxdrift][secondbeforedrift+maxdrift][thirdbeforedrift+maxdrift]+=upmessage[i][xi] *\
                                        ZProbTable[firstbeforedrift+maxdrift][firstdrift-firstbeforedrift+1]\
                                        [secondbeforedrift+maxdrift][seconddrift-secondbeforedrift+1]\
                                        [thirdbeforedrift+maxdrift][thirddrift-thirdbeforedrift+1][nofs]*\
                                        beforeleftmessage

        #normalization process
        if np.sum(leftmessage[i-1])==0:
            print(leftmessage[i-1],i)
        leftmessage[i-1]=sumone(leftmessage[i-1])
        #if i<5:
            #print(rightmessage[i][maxdrift][maxdrift])
        if i<=len(drift[0]) and i<=len(drift[1]) and i<=len(drift[2]):
            maxprob=argmax_ndim(leftmessage[i-1])
            print(drift[0][i-1],drift[1][i-1],drift[2][i-1],maxprob)

    #print(leftmessage)
    return leftmessage

def SetUprightmessage(int block_length,list drift,int maxdrift,int nofc,int period,double pd,double pi,double ps,cnp.ndarray upmessage,list Z, cnp.ndarray ZProbTable):
    cdef int i,limit,xi,firstdrift,seconddrift,thirddrift,firstnextdrift,secondnextdrift,thirdnextdrift,messagesize,nofs,Zindex
    cdef double beforerightmessage
    cdef cnp.ndarray rightmessage
    messagesize=2*maxdrift+1
    rightmessage=np.zeros((block_length,messagesize,messagesize,messagesize),dtype=float)
    rightmessage[0][maxdrift][maxdrift][maxdrift] = 1
    #print("deciding rightmessage")
    for i in range(0, block_length - 1):
        if i < maxdrift:
            limit=i
            nextlimit=i+1
        else:
            limit=maxdrift
            nextlimit=maxdrift
        for firstdrift in range(-limit, limit+1):
            if i - firstdrift > len(Z[0]) - 1:
                continue;
            for firstnextdrift in range(firstdrift - 1, firstdrift + 2):
                if i - firstnextdrift > len(Z[0]) - 1 or firstnextdrift > nextlimit or firstnextdrift <-nextlimit:
                    continue;
                for seconddrift in range(-limit, limit+1):
                    if i - seconddrift > len(Z[1]) - 1:
                        continue;
                    for secondnextdrift in range(seconddrift - 1, seconddrift + 2):
                        if i - secondnextdrift > len(Z[1]) - 1 or secondnextdrift > nextlimit or secondnextdrift <-nextlimit:
                            continue;
                        for thirddrift in range(-limit, limit+1):
                            if i - thirddrift > len(Z[2]) - 1:
                                continue;
                            for thirdnextdrift in range(thirddrift - 1, thirddrift + 2):
                                if i - thirdnextdrift > len(Z[2]) - 1 or thirdnextdrift > nextlimit or thirdnextdrift <-nextlimit:
                                    continue;
                                beforerightmessage= rightmessage[i][firstdrift + maxdrift][seconddrift +maxdrift][thirddrift+maxdrift]
                                #print(thirdnextdrift+maxdrift,rightmessage[i+1][firstnextdrift+maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift])
                                nowdrift=np.array([firstdrift,seconddrift,thirddrift])
                                nextdrift=np.array([firstnextdrift,secondnextdrift,thirdnextdrift])
                                for xi in range(4):
                                    nofs=0
                                    for Zindex in range(3):
                                        for possibles in range(nowdrift[Zindex]-nextdrift[Zindex]+1):
                                            if i-nowdrift[Zindex]+possibles>len(Z[Zindex])-1:
                                                continue;
                                            if str(xi)!=Z[Zindex][i-nowdrift[Zindex]+possibles]:
                                                nofs=+1
                                    

                                    rightmessage[i + 1][firstnextdrift + maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]+=upmessage[i][xi] *\
                                    ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1]\
                                    [seconddrift+maxdrift][secondnextdrift-seconddrift+1]\
                                    [thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*\
                                    beforerightmessage

        #normalization processi
        rightmessage[i+1]=sumone(rightmessage[i+1])
        #if i%10==1:
            #print(rightmessage[i][maxdrift][maxdrift])
        #maxprob=argmax_ndim(rightmessage[i+1])
        #print(drift[0][i+1],drift[1][i+1],drift[2][i+1],maxprob)
    #print(rightmessage)
    return rightmessage

def SetUpDownmessage(int block_length,int maxdrift,int period,double pd,double pi,double ps,list Z,cnp.ndarray leftmessage,cnp.ndarray rightmessage,cnp.ndarray ZProbTable):
    cdef int i,limit,xi,firstdrift,firstnextdrift,seconddrift,secondnextdrift,thirddrift,thirdnextdrift,nofs,possibles,Zindex
    #cdef list down=[[0 for i in range(4)] for j in range(block_length)]
    cdef cnp.ndarray downmessage
    downmessage=np.zeros((block_length,4),dtype=float)
    print('deciding downmessage')
    for i in range(block_length):
        if i<maxdrift:
            limit=i
            nextlimit=i+1
        else:
            limit=maxdrift
            nextlimit=maxdrift
        for firstdrift in range(-limit,limit+1):
            if i-firstdrift>len(Z[0])-1:
                continue;
            for firstnextdrift in range(firstdrift-1,firstdrift+2):
                if i-firstnextdrift>len(Z[0])-1 or firstnextdrift>nextlimit or firstnextdrift<-nextlimit:
                    continue;
                for seconddrift in range(-limit,limit+1):
                    if i-seconddrift>len(Z[1])-1:
                        continue;
                    for secondnextdrift in range(seconddrift-1,seconddrift+2):
                        if i-secondnextdrift>len(Z[1])-1 or secondnextdrift>nextlimit or secondnextdrift<-nextlimit:
                            continue;
                        for thirddrift in range(-limit,limit+1):
                            if i-thirddrift>len(Z[2])-1:
                                continue;
                            for thirdnextdrift in range(thirddrift-1,thirddrift+2):
                                if i-thirdnextdrift>len(Z[2])-1 or thirdnextdrift>nextlimit or thirdnextdrift<-nextlimit:
                                    continue;
                                nowdrift=np.array([firstdrift,seconddrift,thirddrift])
                                nextdrift=np.array([firstnextdrift,secondnextdrift,thirdnextdrift])
                                for xi in range(4):
                                    nofs=0
                                    for Zindex in range(3):
                                        for possibles in range(nowdrift[Zindex]-nextdrift[Zindex]+1):
                                            if i-nowdrift[Zindex]+possibles>len(Z[Zindex])-1:
                                                continue;
                                            if str(xi)!=Z[Zindex][i-nowdrift[Zindex]+possibles]:
                                                nofs+=1 
                                    #print(nofs,downmessage[i][xi])
                                    downmessage[i][xi]+=ZProbTable[firstdrift+maxdrift][firstnextdrift-firstdrift+1]\
                                                        [seconddrift+maxdrift][secondnextdrift-seconddrift+1]\
                                                        [thirddrift+maxdrift][thirdnextdrift-thirddrift+1][nofs]*\
                                                        rightmessage[i][firstdrift+maxdrift][seconddrift+maxdrift][thirddrift+maxdrift]*\
                                                        leftmessage[i][firstnextdrift+maxdrift][secondnextdrift+maxdrift][thirdnextdrift+maxdrift]
        downmessage[i] = sumone(downmessage[i])
    #print(downmessage)
    return downmessage

def sumone(x,axis=None):
    a=np.sum(x)
    result=x/a
    return result

def argmax_ndim(arg_array):
    return np.unravel_index(arg_array.argmax(), arg_array.shape)

