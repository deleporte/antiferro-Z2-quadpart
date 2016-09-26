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
    