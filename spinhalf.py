import scipy.sparse
import numpy as np

def save_sparse_csr(filename,array):
    np.savez(filename, data=array.data, indices=array.indices,
             indptr=array.indptr, shape=array.shape)

def load_sparse_csr(filename):
    loader=np.load(filename)
    return csr_matrix((loader['data'],loader['indices'],loader['indptr']),
                       loader['shape'])

def save_sparse_coo(filename,array):
    np.savez(filename, data=array.data, row=array.row,
             col=array.col,shape=array.shape)

def load_sparse_coo(filename):
    loader=np.load(filename)
    return coo_matrix((loader['data'],(loader['row'],loader['col'])),
                       loader['shape'])

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
                    M[i,number(2,statemod)]=1.
    return M

#def dimermatrix(delta=1.):
    
def storedata():
    deltas=0.1*np.array(range(21))-1.
    for delta in deltas:
        filename=matrix+str(delta)+'.npz'
        M=fullmatrix(delta)
        save_sparse_coo(filename,M)
    
