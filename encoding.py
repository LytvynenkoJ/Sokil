from ntt import *
from sokil_constants import *

def encode_i8(f, w):
    l=w*len(f)
    if l%8!=0:
        l=l+(8-l%8)
    byte=[0 for i in range(l)]
    for i in range(len(f)):
        src = f[i]
        if src<0:
            byte[w*i] = 1
            src += 2**(w-1)
        for j in range(w-1):
            byte[w*(i+1)-j - 1] = src%2
            src//=2
    j = 0
    result=""
    while j < len(byte)//8:
        src=0
        j+=1
        for i in range(8):
            src+=(2**i)*byte[8*j-i-1]
        result+=str(hex(src))
        result+=" "
    return j, result

def decode_i8(r, w):
    #r = f.split()
    for i in range(len(r)):
        r[i] = int(bin(int(r[i], base=16)), base=2)
    rb = [0 for i in range(8*len(r))]
    for i in range(len(r)):
        src = r[i]
        for j in range(8):
            rb[(i+1)*8 - j - 1] = src%2
            src//=2
    result=[0 for i in range(len(rb)//w)]
    for i in range(len(rb)//w):
        src = 0
        for j in range(w-1):
            src+=(2**j)*rb[w*(i+1)-j-1]
        if rb[w*i]==1:
            src-=2**(w-1)
        result[i] = src
    return result

def recoveryG(f, g, F):
    Fg = intt(mul_ntt(ntt(F), ntt(g)))
    for i in range(len(Fg)):
        Fg[i]+=q
    success, Fgf = div_ntt(ntt(Fg), ntt(f))
    if success:
        G = intt(Fgf)
        for i in range(len(G)):
            if G[i]>q//2:
                G[i]-=q
    return G

def encode_i16(f, w=LOG2_Q):
    l=w*len(f)
    if l%8!=0:
        l=l+(8-l%8)
    byte=[0 for i in range(l)]
    for i in range(len(f)):
        src = f[i]
        for j in range(w):
            byte[w*(i+1)-j - 1] = src%2
            src//=2
    j = 0
    result=""
    while j < len(byte)//8:
        src=0
        j+=1
        for i in range(8):
            src+=(2**i)*byte[8*j-i-1]
        result+=str(hex(src))
        result+=" "
    return j, result

def decode_i16(r, w=LOG2_Q):
    #r = f.split()
    for i in range(len(r)):
        r[i] = int(bin(int(r[i], base=16)), base=2)
    rb = [0 for i in range(8*len(r))]
    for i in range(len(r)):
        src = r[i]
        for j in range(8):
            rb[(i+1)*8 - j - 1] = src%2
            src//=2
    result=[0 for i in range(len(rb)//w)]
    for i in range(len(rb)//w):
        for j in range(w):
            result[i]+=(2**j)*rb[w*(i+1)-j-1]
    return result