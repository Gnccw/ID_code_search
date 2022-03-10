import math

def EK(N,L,logQ):
    b=N*pow(2,L-2*N+logQ)*(1-pow(2,L-2*N-1))
    return pow(math.e,-b)
LNQ=[[16,27.228,3],[24,40.618,5],[32,48.203,13],[48,72.618,20],[64,96.203,28]]
for i in LNQ:
    print(EK(i[0],i[1],i[2]))


def C(N,L,logQ):
    b=N*pow(2,L-2*N+logQ)
    c=pow(math.e,-1*b)
    ddict={'16':14,'24':15,'32':16,'48':19,'64':22,}
    return (3*N-logQ+math.log2(1-c)-math.log2(ddict[str(N)]/2))

for i in LNQ:
    print(C(i[0],i[1],i[2]))

