from sokil_constants import *

modulo = 65536
Q = 12287

def mod_mul(x,y):
    z=x*y
    w=(z*Q)%modulo #(mod 2^{16})
    w*=q
    d=(z+w)//modulo
    d%=q
    return d

def NTT(a,n):
    t = n
    m = 1
    while m < n:
        ht = t // 2
        u = 0
        j1 = 0
        while u < m:
            s = NTT_tabl[m + u]
            j2 = j1 + ht
            j = j1
            while j < j2:
                w = a[j]
                v = mod_mul(a[j + ht], s)
                a[j] = (w + v)%q
                a[j + ht] = (w - v)%q
                j = j + 1
            u+=1
            j1=j1 + t
        m = m * 2
        t = ht
    return a

def InvNTT(a,n):
    t = 1
    m = n
    dt=1
    while m > 1:
        hm = m // 2
        dt *= 2
        r = 0
        j1 = 0
        while r < hm:
            j2 = j1 + t
            s = NTT_tabl_[hm + r]
            j = j1
            while j < j2:
                u = a[j]
                v = a[j + t]
                a[j] = (u + v) % q
                w = (u - v) % q
                a[j + t] = mod_mul(w, s)
                j = j + 1
            j1 = j1 + dt
            r = r + 1
        t = dt
        m = hm
    for m in range(n):
        a[m] = mod_mul(a[m], NI[n])
    return a

def mod_div(x, y):
    R2=10952 #2^{32} mod q
    y0 = mod_mul(y, R2)
    y = y0
    y2=y
    for i in range(2):
        for i in range(4):
            y = mod_mul(y,y)
        y = mod_mul(y,y0)
    y0 = y              #^273
    y = mod_mul(y,y)    #^576
    y = mod_mul(y,y)    #^1092
    y = mod_mul(y,y0)   #^1365
    y = mod_mul(y,y)    #^2730
    y = mod_mul(y,y0)   #^3003
    y = mod_mul(y,y)    #^6006
    y = mod_mul(y,y)    #^12012
    y = mod_mul(y, y0)  #^12285
    y2 = mod_mul(y2,y2)         #y^{2}
    y = mod_mul(y2,y)   #^12287=^{-1}
    z = mod_mul (y, x)
    return z

def poly_montydiv_ntt(a,b,n):
    c=[0 for i in range(n)]
    for i in range(n):
        if b[i]==0:
            return False, c
        c[i] = mod_div(a[i], b[i])
    return True, c

def div_ntt(f_ntt, g_ntt):
    assert len(f_ntt) == len(g_ntt)
    deg = len(f_ntt)
    if any(elt == 0 for elt in g_ntt):
        return False, []
    return True, [(f_ntt[i] * inv_mod_q[g_ntt[i]]) % q for i in range(deg)]