import scipy.linalg
import scipy.sparse
import numpy

def buildSymbolMatrix(Nx=4,Ny=4):
    if Ny==1:
        graph = buildChain(Nx)
    else:
        graph = buildGrid(Nx,Ny)
    Adj = 0.5*networkx.adjacency_matrix(graph)
    if Ny==1:
        M = Adj+ scipy.sparse.identity(Nx*Ny)
        N = -Adj+ scipy.sparse.identity(Nx*Ny)
    else:
        M = Adj+ 2.*scipy.sparse.identity(Nx*Ny)
        N = -Adj+ 2.*scipy.sparse.identity(Nx*Ny)
        
    Z = scipy.sparse.identity(Nx*Ny)
    Z = Z - Z
    return scipy.sparse.vstack([scipy.sparse.hstack([M,Z]),
                                scipy.sparse.hstack([Z,N])])

def buildJSymbolMatrix(Nx=4,Ny=4):
    if Ny==1:
        graph = buildChain(Nx)
    else:
        graph = buildGrid(Nx,Ny)
    Adj = 0.5*networkx.adjacency_matrix(graph)
    if Ny==1:
        M = Adj+ 2.*scipy.sparse.identity(Nx*Ny)
        N = Adj- 2.*scipy.sparse.identity(Nx*Ny)
    else:
        M = Adj+ 2.*scipy.sparse.identity(Nx*Ny)
        N = Adj- 2.*scipy.sparse.identity(Nx*Ny)
    Z = scipy.sparse.identity(Nx*Ny)
    Z = Z - Z
    return scipy.sparse.vstack([scipy.sparse.hstack([Z,M]),
                                scipy.sparse.hstack([N,Z])])

def charac(M):
    N=M.todense()
    values = numpy.imag(scipy.linalg.eigvals(N))
    mu=0
    for v in values:
        if v>0:
            mu+=v
    return mu
    
    
    
    
