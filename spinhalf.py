import scipy.sparse
import numpy as np
import scipy.sparse.linalg

def save_sparse_csr(filename,array):
    np.savez(filename, data=array.data, indices=array.indices,
             indptr=array.indptr, shape=array.shape)

def load_sparse_csr(filename):
    loader=np.load(filename)
    return scipy.sparse.csr_matrix((loader['data'],loader['indices'],
                                    loader['indptr']),loader['shape'])

def save_sparse_coo(filename,array):
    np.savez(filename, data=array.data, row=array.row,
             col=array.col,shape=array.shape)

def load_sparse_coo(filename):
    loader=np.load(filename)
    return coo_matrix((loader['data'],(loader['row'],loader['col'])),
                       loader['shape'])

def mean(j,vector):
    #returns <Sxj>,<Szj>
    mn=[0,0]
    for i in range(2**16):
        state=digits(2,i)
        while len(state)<16:
            state.append(0)
        varstate=state[:]
        varstate[j]=1-varstate[j]
        mn[0]+=vector[i]*vector[number(2,varstate)]
        #mean[1]-=vector[i]*vector[number(2,varstate)] #?
        mn[1]+=vector[i]*vector[i]*(2*state[j]-1)
    return mn

def correlations(j1,j2,vector):
    #by definition it is <SiSj>-<Si><Sj>
    #it is not a value, but a 2x2 matrix
    correl=np.array([[-mean(j1,vector)[0]*mean(j2,vector)[0],
                      -mean(j1,vector)[0]*mean(j2,vector)[1]],
                     [-mean(j1,vector)[1]*mean(j2,vector)[0],
                      -mean(j1,vector)[1]*mean(j2,vector)[1]]])
    for i in range(2**16):
        state=digits(2,i)
        while len(state)<16:
            state.append(0)
        varstate1=state[:]
        varstate1[j1]=1-varstate1[j1]
        varstate2=state[:]
        varstate2[j2]=1-varstate2[j2]
        varstate3=state[:]
        varstate3[j1]=1-varstate3[j1]
        varstate3[j2]=1-varstate3[j2]
        correl[0,0]+=vector[i]*vector[number(2,varstate3)]
        correl[0,1]+=vector[i]*vector[number(2,varstate1)]*(2*state[j2]-1)
        correl[1,0]+=vector[i]*vector[number(2,varstate2)]*(2*state[j1]-1)
        correl[1,1]+=vector[i]*vector[i]*(2*state[j1]-1)*(2*state[j2]-1)
    return correl
        

def density(j1,j2,vector):
    d=0.
    for i in range(2**14):
        dig = digits(2,i)
        while len(dig)<14:
            dig.append(0)
        state=dig[:j1]
        state.append(1)
        state=state+dig[j1:j2]
        state.append(1)
        state=state+dig[j2:]
        d += vector[number(2,state)]**2

        state[j1]=0
        state[j2]=0
        d += vector[number(2,state)]**2

        state[j1]=0
        state[j2]=1
        d -= vector[number(2,state)]**2
        
        state[j1]=1
        state[j2]=0
        d -= vector[number(2,state)]**2
        
    return d

def number(base,digits):
    result=0
    temp=digits[:]
    while len(temp) != 0:
        result *= base
        result += temp.pop()
    return result

def digits(base,number):
    digits=[]
    N=int(number)
    while N > 0:
        digits.append(N-base*int(N/base))
        N /= base
    return digits

def fullmatrix(delta=1.):
    M = 0.5*scipy.sparse.identity(65536)
    M = M-M
    for i in range(65536):
        print i
        state = digits(2,i)
        while len(state)<16:
            state.append(0)
        nb=range(4)
        for j in range(16):
            nb[0]=j+1
            if not nb[0]%4:
                nb[0]-=4
            nb[1]=j-1
            if not j%4:
                nb[1]+=4
            nb[2]=j+4
            if nb[2]>=16:
                nb[2]-=16
            nb[3]=j-4
            if nb[3]<0:
                nb[3]+=16
            for k in range(4):
                if state[j]==state[nb[k]]:
                    M[i,i]-=delta
                else:
                    M[i,i]+=delta
                    statemod=state[:]
                    statemod[j]=1-statemod[j]
                    statemod[nb[k]]=1-statemod[nb[k]]
                    M[i,number(2,statemod)]=-1.
    return M

#def dimermatrix(delta=1.):

def transformmatrix(M,delta):
    for i in range(65536):
        M[i,i]=-1.*float(delta)*M[i,i]
    
def storevectors(delta,Nb_values=100):
    M=load_sparse_csr('matrix-1.0.npz')
    transformmatrix(M,delta)
    val,vect=scipy.sparse.linalg.eigsh(M,Nb_values)
    filename='vectors'+str(delta)+'.npz'
    np.savez(filename,values=val,vectors=vect)

def loadvectors(delta):
    filename='vectors'+str(delta)+'.npz'
    loader=np.load(filename)
    return loader['values'],loader['vectors']
    
def storedata():
    deltas=0.1*np.array(range(21))-1.
    for delta in deltas:
        filename='matrix'+str(delta)+'.npz'
        M=fullmatrix(delta)
        save_sparse_csr(filename,M)
    
