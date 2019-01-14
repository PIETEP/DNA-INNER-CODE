def setprobtable(maxdrift,pd,pi,ps):


    DriftProbTable=setdriftprobtable(maxdrift,pd,pi)
    JointDriftProbTable=setjointdriftprobtable(DriftProbTable,maxdrift) 
    
    ZProbTable=setZprobtable(JointDriftProbTable,maxdrift,ps)
    return ZProbTable

def setdriftprobtable(maxdrift,pd,pi):
    DriftProbTable=[[0 for j in range(3)] for i in range(2*maxdrift+1)]
    tempprob=0
    for drift in range(-maxdrift,maxdrift+1):
        for nextdrift in range(3):
            if nextdrift==0:
                if drift==-maxdrift:
                    continue;
                else:
                    tempprob=pi
            if nextdrift==1:
                tempprob=1-pd-pi
            else:
                if drift==maxdrift:
                    continue;
                else:
                    tempprob=pd
            DriftProbTable[drift+maxdrift][nextdrift]=tempprob
    return DriftProbTable

def setjointdriftprobtable(DriftProbTable,maxdrift):
    JointDriftProbTable=[[[[[[0 for thirdnextdrift in range(3)] for thirddrift in range(2*maxdrift+1)] for secondnextdrift in range(3)] for seconddrift in range(2*maxdrift+1)\
                        ] for firstnextdrift in range(3)] for firstdrift in range(2*maxdrift+1)]
    pfirst=0
    psecond=0
    pthird=0
    for firstdrift in range(2*maxdrift+1):
        for firstnextdrift in range(3):
            if firstdrift+firstnextdrift==0 or firstdrift+firstnextdrift==2*maxdrift+2:
                continue;
            pfirst=DriftProbTable[firstdrift][firstnextdrift]
            for seconddrift in range(2*maxdrift+1):
                for secondnextdrift in range(3):
                    if seconddrift+secondnextdrift==0 or seconddrift+secondnextdrift==2*maxdrift+2:
                        continue;
                    psecond=DriftProbTable[seconddrift][secondnextdrift]*pfirst
                    for thirddrift in range(2*maxdrift+1):
                        for thirdnextdrift in range(3):
                            if thirddrift+thirdnextdrift==0 or thirddrift+thirdnextdrift==2*maxdrift+2:
                                continue;
                            pthird=DriftProbTable[thirddrift][thirdnextdrift]
                            JointDriftProbTable[firstdrift][firstnextdrift][seconddrift][secondnextdrift][thirddrift][thirdnextdrift]=psecond*pthird

    return JointDriftProbTable


def setZprobtable(JointDriftProbTable,maxdrift,ps):
    ZProbTable=[[[[[[[0 for nofs in range(7)] for thirdnextdrift in range(3)] for thirddrift in range(2*maxdrift+1)] for secondnextdrift in range(3)] for seconddrift in range(2*maxdrift+1)\
                        ] for firstnextdrift in range(3)] for firstdrift in range(2*maxdrift+1)] 
    for firstdrift in range(2*maxdrift+1):
        for firstnextdrift in range(3):
            if firstdrift+firstnextdrift==0 or firstdrift+firstnextdrift==2*maxdrift+2:
                continue;
            for seconddrift in range(2*maxdrift+1):
                for secondnextdrift in range(3):
                    if seconddrift+secondnextdrift==0 or seconddrift+secondnextdrift==2*maxdrift+2:
                        continue;
                    for thirddrift in range(2*maxdrift+1):
                        for thirdnextdrift in range(3):
                            if thirddrift+thirdnextdrift==0 or thirddrift+thirdnextdrift==2*maxdrift+2:
                                continue;
                            delinscount=[0 for count in range(3)]
                            delinscount[firstnextdrift]+=1
                            delinscount[secondnextdrift]+=1
                            delinscount[thirdnextdrift]+=1
                            DriftProb=JointDriftProbTable[firstdrift][firstnextdrift][seconddrift][secondnextdrift][thirddrift][thirdnextdrift]
                            for nofs in range(7):
                                if nofs>2*delinscount[0]+delinscount[1]:
                                    Zprob=0
                                else:
                                    Zprob=pow(1-ps,2*delinscount[0]+delinscount[1]-nofs)*pow(ps/3,nofs)
                                ZProbTable[firstdrift][firstnextdrift][seconddrift][secondnextdrift][thirddrift][thirdnextdrift][nofs]=DriftProb*Zprob
    return ZProbTable
