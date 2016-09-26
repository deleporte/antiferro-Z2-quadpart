import numpy as np
import scipy.sparse

def number(base,digits):
    result=0
    while len(digits) != 0:
        result *= base
        result += digits.pop()
    return result

def digits(base,number):
    digits=[]
    N=int(number)
    while N > 0:
        digits.append(N-base*int(N/base))
        N /= base
    return digits[::-1]

def buildOpMatrix(Nx=4,Ny=4,cutoff=2):
    M=0.5*scipy.sparse.identity(cutoff**(Nx*Ny))
    M=M-M
    for i in range(cutoff**(Nx*Ny)):
        state = digits(cutoff,i)
        #diagonal term
        diag = 0
        for j in range(Nx*Ny):
            diag += 2.*(state[j]+1)
        M[i,i]=diag
        #off-diagonal terms
        nb=[0,0,0,0]
        for j in range(Nx*Ny):
            #who are the neighbours of j?
            nb[0]=j+1
            if not nba%Nx:
                nb[0] -= Nx
            nb[1]=j-1
            if not j%Nx:
                nb[1] += Nx
            nb[2]=j+Nx
            if nb[2] >= Nx*Ny:
                nb[2] -= Nx*Ny
            nb[3]=j-Nx
            if nb[3] <0:
                nb[3] += Nx*Ny
            otherstate=state[:]
            for k in range(4):
                otherstate[j]+=1
                otherstate[nb[k]]+=1
                value = 0.5*np.sqrt(otherstate[nb[k]]*otherstate[j])
                M[i,number(cutoff,otherstate)]=value
                otherstate[j]-=1
                otherstate[nb[k]]-=1
                value = 0.5*np.sqrt(otherstate[nb[k]]*otherstate[j])
                otherstate[j]-=1
                otherstate[nb[k]]-=1
                M[i,number(cutoff,otherstate)]=value
                otherstate[j]+=1
                otherstate[nb[k]]+=1
    return M
                
