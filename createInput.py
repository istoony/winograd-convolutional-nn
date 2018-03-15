
# coding: utf-8

# In[63]:


import random
import sys


# In[64]:


"""
create fname (default input.txt) of type int (default) or float of w*h*ch size
w = width
h = height    --only square is tested
ch = channels
"""

def createIMG(w, h, ch, t='i', fname="input.txt"):
    if t=='f':
        with open(fname, "w") as f:
            for i in range(0, w*h*ch):
               f.write(str(random.random())+" ")
    else:
        with open(fname, "w") as f:
            for i in range(0, w*h*ch):
               f.write(str(random.randint(0,10))+" ")
        
              


# In[65]:


w = int(sys.argv[1])
h = int(sys.argv[2])
ch = int(sys.argv[3])

if len(sys.argv)<5:
    createIMG(w, h, ch)

if len(sys.argv)==5:
    if len(sys.argv[4])>1:
        f = sys.argv[4]
        createIMG(w, h, ch, fname=f)
    else:
        ty = sys.argv[4]
        createIMG(w, h, ch, t=ty)
    
if len(sys.argv)==6:  
    if len(sys.argv[4])>1:
        f = sys.argv[4]
        ty = sys.argv[5]
    else:
        ty = sys.argv[4]
        f = sys.argv[5]
    createIMG(w, h, ch, t=ty, fname=f)

