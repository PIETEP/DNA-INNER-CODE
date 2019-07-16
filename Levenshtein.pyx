def levenshtein_distance(a, b):
    #distance,deletion,insertion,substitution
    m = [ [[0,0,0,0] for j in range (len(b) + 1) ]for i in range(len(a) + 1) ]
    #print(m[0][0])
    for i in range(len(a) + 1):
        m[i][0][0] = i
        m[i][0][1] = i

    for j in range(len(b) + 1):
        m[0][j][0] = j
        m[0][j][2] = j

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                x = 0
            else:
                x = 1
            nextm=min(m[i - 1][j][0] + 1, m[i][ j - 1][0] + 1, m[i - 1][j - 1][0] + x)
            m[i][j][0]=nextm
            if nextm == m[i-1][j-1][0]+x:
                m[i][j][1]=m[i-1][j-1][1]
                m[i][j][2]=m[i-1][j-1][2]
                m[i][j][3]=m[i-1][j-1][3]
                if x==1:
                    m[i][j][3]+=1
            if nextm ==m[i-1][j][0]+1:
                if nextm != m[i-1][j-1][0]+x:
                    m[i][j][1]=m[i-1][j][1]+1
                    m[i][j][2]=m[i-1][j][2]
                    m[i][j][3]=m[i-1][j][3]
                else:
                    if m[i][j][1]+m[i][j][2]>m[i-1][j][1]+m[i-1][j][2]+1:
                        m[i][j][1]=m[i-1][j][1]+1
                        m[i][j][2]=m[i-1][j][2]
                        m[i][j][3]=m[i-1][j][3]
            if nextm == m[i][j-1][0]+1:
                if nextm !=m[i-1][j-1][0]+x and nextm!=m[i-1][j][0]+1:
                    m[i][j][1]=m[i][j-1][1]
                    m[i][j][2]=m[i][j-1][2]+1
                    m[i][j][3]=m[i][j-1][3]
                else:
                    if m[i][j][1]+m[i][j][2]>m[i][j-1][1]+m[i][j-1][2]+1:
                        m[i][j][1]=m[i][j-1][1]
                        m[i][j][2]=m[i][j-1][2]+1
                        m[i][j][3]=m[i][j-1][3]

    # print m
    #nosub=int(m[-1][-1])
    #print(m[-1][-1],m[-1][-1]-nosub,round(1-round((m[-1][-1]-nosub),3),3))
    #Levenshtein=((round(1-round((m[-1][-1]-nosub),3),3))*1000)%1000+nosub
    return m[-1][-1][0],m[-1][-1][2],m[-1][-1][1],m[-1][-1][3]
