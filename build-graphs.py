import networkx

def buildChain(N=4):
    graph=networkx.Graph()
    for k in range(N-1):
        graph.add_edge(k,k+1)
    graph.add_edge(N-1,0)
    return graph

def buildGrid(Nx=4,Ny=4):
    graph=networkx.Graph()
    for k in range(Nx-1):
        for l in range(Ny-1):
            graph.add_edge((k,l),(k+1,l))
            graph.add_edge((k,l),(k,l+1))
        graph.add_edge((k,Ny-1),(k,0))
        graph.add_edge((k,Ny-1),(k+1,Ny-1))
    for l in range(Ny-1):
        graph.add_edge((Nx-1,l),(0,l))
        graph.add_edge((Nx-1,l),(Nx-1,l+1))
    graph.add_edge((Nx-1,Ny-1),(0,Ny-1))
    graph.add_edge((Nx-1,Ny-1),(Nx-1,0))
    return graph

            
