
# coding: utf-8

# In[1]:


import random
import sys


# In[3]:


"""
create fname (default kernels.txt) of type int (default) or float of w*h*ch*n size
w = width
h = height    --only square is tested
ch = channels
n = number of kernels
"""
def createKernels(w, h, ch, n, t='i', fname="kernels.txt"):
    if t=='f':
        with open(fname, "w") as f:
            for i in range(0, w*h*ch*n):
               f.write(str(10*random.random())+" ")
    else:
        with open(fname, "w") as f:
            for i in range(0, w*h*ch*n):
               f.write(str(random.randint(0,10))+" ")
        


# In[4]:


w = int(sys.argv[1])
h = int(sys.argv[2])
ch = int(sys.argv[3])
n = int(sys.argv[4])

if len(sys.argv)<6:
    createKernels(w, h, ch, n)

if len(sys.argv)==6:
    if len(sys.argv[5])>1:
        f = sys.argv[5]
        createKernels(w, h, ch, n, fname=f)
    else:
        ty = sys.argv[5]
        createKernels(w, h, ch, n, t=ty)
    
if len(sys.argv)==7:  
    if len(sys.argv[5])>1:
        f = sys.argv[5]
        ty = sys.argv[6]
    else:
        ty = sys.argv[5]
        f = sys.argv[6]
    createKernels(w, h, ch, n, t=ty, fname=f)

