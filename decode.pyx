import factor_graph
import joint_graph
import Levenshtein
import numpy as np
cimport numpy as cnp

def decode(block_length,CopyZ,cordword,decordtimes,maxdrift,nofcopy,NextPrior,period,pid,ps,priorlist):
    for decode_index in range(decordtimes):
        print('decode[%d/%d]' %(decode_index+1,decordtimes))
        for copy_index in range(nofcopy):
            #print(CopyZ[copy_index])
            normalized_out_probability=factor_graph.factor_graph(block_length,CopyZ[copy_index],maxdrift,NextPrior,period,pid,ps,priorlist)
            if copy_index!=0:
                out_probability_list=np.dstack([out_probability_list,normalized_out_probability])
            else:
                out_probability_list=normalized_out_probability

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


    flag=0
    decodeerrorcount=0
    decode_word_error_count=0
    for i in range(NextPrior.shape[0]):
        predicted_base=argmax_ndim(NextPrior[i])
        #print(predicted_base[0])
        decode_word.append(bin(predicted_base[0]))
        if decode_word[i]!=cordword[i]:
            decodeerrorcount+=1
            if flag==0:
                flag=1
                decode_word_error_count+=1

    return decodeerrorcount,decode_word_error_count

def joint_decode(block_length,cordword,Copydrift,CopyZ,maxdrift,nofcopy,NextPrior,period,pd,pi,ps,ZProbTable):
    
    normalized_out_probability=joint_graph.factor_graph(block_length,Copydrift,CopyZ,maxdrift,nofcopy,NextPrior,period,pd,pi,ps,ZProbTable)
    
    decode_word=""
    correct_probability=[]
    flag=0
    decodeerrorcount=0
    decode_word_error_count=0
    for i in range(normalized_out_probability.shape[0]):
        predicted_base=argmax_ndim(normalized_out_probability[i])
        #print(predicted_base[0])
        decode_word+=str(predicted_base[0])
        if decode_word[i]!=cordword[i]:
            decodeerrorcount+=1
            print(normalized_out_probability[i][predicted_base],normalized_out_probability[i][int(cordword[i])],i)
            if flag==0:
                flag=1
                decode_word_error_count+=1
        elif flag==0:
            correct_probability.append(normalized_out_probability[i][predicted_base])
    if flag==0:
        correct_probability.sort()
        correctlist=correct_probability[0:10]
    else:
        correctlist=[0]*10
    correctlist=np.array(correctlist)

    leven,ins,dele,sub=Levenshtein.levenshtein_distance(cordword,decode_word)
    return decodeerrorcount,decode_word_error_count,correctlist,leven,ins,dele,sub
    
def argmax_ndim(arg_array):
    return np.unravel_index(arg_array.argmax(), arg_array.shape)


