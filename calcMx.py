
# coding: utf-8

# In[15]:


from sympy import symbols, Matrix, Poly, zeros, eye, Indexed, simplify, IndexedBase, init_printing, pprint
from operator import mul
from functools import reduce
import numpy as np
import sys


def At(a,m,n):
    return Matrix(m, n, lambda i,j: a[i]**j)

def A(a,m,n):
    return At(a, m-1, n).row_insert(m-1, Matrix(1, n, lambda i,j: 1 if j==n-1 else 0))

def T(a,n):
    return Matrix(Matrix.eye(n).col_insert(n, Matrix(n, 1, lambda i,j: -a[i]**n)))

def Lx(a,n):
    x=symbols('x')
    return Matrix(n, 1, lambda i,j: Poly((reduce(mul, ((x-a[k] if k!=i else 1) for k in range(0,n)), 1)).expand(basic=True), x))

def F(a,n):
    return Matrix(n, 1, lambda i,j: reduce(mul, ((a[i]-a[k] if k!=i else 1) for k in range(0,n)), 1))

def Fdiag(a,n):
    f=F(a,n)
    return Matrix(n, n, lambda i,j: (f[i,0] if i==j else 0))

def FdiagPlus1(a,n):
    f = Fdiag(a,n-1)
    f = f.col_insert(n-1, zeros(n-1,1))
    f = f.row_insert(n-1, Matrix(1,n, lambda i,j: (1 if j==n-1 else 0)))
    return f

def L(a,n):
    lx = Lx(a,n)
    f = F(a, n)
    return Matrix(n, n, lambda i,j: lx[i, 0].nth(j)/f[i]).T

def Bt(a,n):
    return L(a,n)*T(a,n)

def B(a,n):
    return Bt(a,n-1).row_insert(n-1, Matrix(1, n, lambda i,j: 1 if j==n-1 else 0))

FractionsInG=0
FractionsInA=1
FractionsInB=2
FractionsInF=3

def cookToomFilter(a,n,r,fractionsIn=FractionsInG):
    alpha = n+r-1
    f = FdiagPlus1(a,alpha)
    if f[0,0] < 0:
        f[0,:] *= -1
    if fractionsIn == FractionsInG:
        AT = A(a,alpha,n).T
        G = (A(a,alpha,r).T/f).T
        BT = f * B(a,alpha).T
    elif fractionsIn == FractionsInA:
        BT = f * B(a,alpha).T
        G = A(a,alpha,r)
        AT = (A(a,alpha,n)).T/f
    elif fractionsIn == FractionsInB:
        AT = A(a,alpha,n).T
        G = A(a,alpha,r)
        BT = B(a,alpha).T
    else:
        AT = A(a,alpha,n).T
        G = A(a,alpha,r)
        BT = f * B(a,alpha).T
    return (AT,G,BT,f)


def filterVerify(n, r, AT, G, BT):

    alpha = n+r-1

    di = IndexedBase('d')
    gi = IndexedBase('g')
    d = Matrix(alpha, 1, lambda i,j: di[i])
    g = Matrix(r, 1, lambda i,j: gi[i])

    V = BT*d
    U = G*g
    M = U.multiply_elementwise(V)
    Y = simplify(AT*M)

    return Y

def convolutionVerify(n, r, B, G, A):

    di = IndexedBase('d')
    gi = IndexedBase('g')

    d = Matrix(n, 1, lambda i,j: di[i])
    g = Matrix(r, 1, lambda i,j: gi[i])

    V = A*d
    U = G*g
    M = U.multiply_elementwise(V)
    Y = simplify(B*M)

    return Y

def calcFNames(n, r):
    at = "AT"+str(n)+"x"+str(r)+".txt"
    g = "G"+str(n)+"x"+str(r)+".txt"
    bt = "BT"+str(n)+"x"+str(r)+".txt"
    return at, g, bt

def showCookToomFilter(a,n,r,fractionsIn=FractionsInG):

    AT,G,BT,f = cookToomFilter(a,n,r,fractionsIn)
    
    fnameAT, fnameG, fnameBT = calcFNames(n, r)
    
    with open(fnameAT, "wb") as fa:
        #print("AT = ")
        #pprint(AT)
        np.savetxt(fa, np.array(AT).astype(np.float64), fmt="%.10f")
        #print("")
        
    with open(fnameG, "wb") as fg:
        #print("G = ")
        #pprint(G)
        np.savetxt(fg, np.array(G).astype(np.float64), fmt="%.10f")
        #print("")
        
    with open(fnameBT, "wb") as fb:
        #print("BT = ")
        #pprint(BT)
        np.savetxt(fb, np.array(BT).astype(np.float64), fmt="%.10f")
        #print("")
    

    #if fractionsIn != FractionsInF:
        #print("FIR filter: AT*((G*g)(BT*d)) =")
        #pprint(filterVerify(n,r,AT,G,BT))
        #print("")

    #if fractionsIn == FractionsInF:
        #print("fractions = ")
        #pprint(f)
        #print("")

def showCookToomConvolution(a,n,r,fractionsIn=FractionsInG):

    AT,G,BT,f = cookToomFilter(a,n,r,fractionsIn)

    B = BT.transpose()
    A = AT.transpose()
    
    print("A = ")
    pprint(A)
    print("")
    
        
    
    print("G = ")
    pprint(G)
    print("")

    print("B = ")
    pprint(B)
    print("")

    if fractionsIn != FractionsInF:
        print("Linear Convolution: B*((G*g)(A*d)) =")
        pprint(convolutionVerify(n,r,B,G,A))
        print("")

    if fractionsIn == FractionsInF:
        print("fractions = ")
        pprint(f)
        print("")
        
        
def calcInterpolaionPoints(n, r):
    p = n+r-2
    l = [0]
    return calcIPRec(l, p-1, 1)
    
    
def calcIPRec(v, n, p):
    if n<=0:
        return v
    v.append(p)
    v.append(-1*p)
    calcIPRec(v, n-2, p+1)
    return v

    


# In[16]:


n = int(sys.argv[1])
r = int(sys.argv[2])

ip = calcInterpolaionPoints(n, r)

#F(n,r)

showCookToomFilter(ip , n, r)

