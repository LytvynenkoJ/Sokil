from gmpy2 import mpfr
from copy import deepcopy

from sokil_constants import *
from common import *
from sokil_fft import *
from ntt import ntt, intt, div_ntt
from gauss import *
from fft import *

def gen_h(g,f,n):
    g_ntt = ntt(g)
    f_ntt = ntt(f)
    success, h = div_ntt(g_ntt, f_ntt)
    if success:
        h = intt(h)
    return success, h

def gen_fgh(n):
    success=False
    while success==False:
        f = poly_mkgauss(n)
        g = poly_mkgauss(n)

        f_copy = deepcopy(f)
        g_copy = deepcopy(g)

        normf=sqnorm(f)
        normg=sqnorm(g)
        if normf+normg>max_norma:
            continue
        
        fft_f = cpoly_fft(f, n); fft_g = cpoly_fft(g, n)
        temp = сpoly_add_muladj_fft(fft_f, fft_g, fft_f, fft_g, n)
        fft_fs = cpoly_adj(fft_f,n); fft_gs = cpoly_adj(fft_g,n)
        if has_Null(temp):
            continue
        for i in range(len(fft_fs)):
            fft_fs[i]*=q
            fft_gs[i]*=q
        f2 = сpoly_ddiv_fft(fft_fs, temp)
        g2 = сpoly_ddiv_fft(fft_gs, temp)
        f=cpoly_inv_fft(f2,n)
        g=cpoly_inv_fft(g2,n)
        
        normf=sqnorm(f)#*(q**2)
        normg=sqnorm(g)#*(q**2)
        if normf+normg>max_norma:
            continue
        
        success, h = gen_h(g_copy, f_copy, n)
    return f_copy, g_copy, h

def gen_half(src, t):
    n2=t//2
    temp_e = [src[2 * i] for i in range(n2)]
    temp_o = [src[2 * i + 1] for i in range(n2)]

    temp_e2 = karamul(temp_e, temp_e)
    temp_o2 = karamul(temp_o, temp_o)
    
    dest=[0]
    for i in range(n2-1):
        dest.append(temp_e2[i + 1] - temp_o2[i])
    dest[0] = temp_e2[0] + temp_o2[n2 - 1]
    return dest

def FG_step1(f,g,n):
    lf=[[mpfr(0) for j in range(n)] for i in range(logn[n]+1)]
    lg=[[mpfr(0) for j in range(n)] for i in range(logn[n]+1)]
    for i in range(n):
        lf[0][i]=f[i]
        lg[0][i]=g[i]
    k = 1
    u = n
    while u >=2:
        lf[k] = gen_half(lf[k - 1], u)
        lg[k] = gen_half(lg[k - 1], u)
        k += 1
        u = u //2
    return lf, lg

def reduce(lF, lG, lf, lg, t):
    sc = max(53, bitsize(min(lf)), bitsize(max(lf)), bitsize(min(lg)), bitsize(max(lg))) - 53

    lf1 = [elt >> (sc) for elt in lf]
    lg1 = [elt >> (sc) for elt in lg]
    
    df1=[]
    dg1=[]
    if t!=2048:
        df1 = fft(lf1)
        dg1 = fft(lg1)
        dtemp1 = add_fft(mul_fft(df1, adj_fft(df1)), mul_fft(dg1, adj_fft(dg1)))
    else:
        df1 = cpoly_fft(lf1, t)
        dg1 = cpoly_fft(lg1, t)
        dtemp1 = сpoly_add_muladj_fft(df1, df1, dg1, dg1, t)
    
    while (1):
        sc1 = max(53, bitsize(min(lF)), bitsize(max(lF)), bitsize(min(lG)), bitsize(max(lG))) - 53
        if sc1 < sc:
            return lF, lG
        
        lF1 = [elt >> (sc1) for elt in lF]
        lG1 = [elt >> (sc1) for elt in lG]
        
        if t!=2048:
            dF1 = fft(lF1)
            dG1 = fft(lG1)
            dtemp2 = add_fft(mul_fft(dF1, adj_fft(df1)), mul_fft(dG1, adj_fft(dg1)))
            dk_fft = div_fft(dtemp2, dtemp1)
            dk = ifft(dk_fft)
        else:
            dF1 = cpoly_fft (lF1, t)
            dG1 = cpoly_fft (lG1, t)
            dtemp2 = сpoly_add_muladj_fft(dF1, df1, dG1, dg1, t)
            dk_fft = сpoly_ddiv_fft(dtemp2, dtemp1)
            dk = cpoly_inv_fft (dk_fft, t)
                
        k = [int(round(dk[i])) for i in range(t)]
        #print("k = " + str(k))
        
        success = True
        for i in range(t):
            if k[i]!=0:
                success=False
        if success:
            return lF, lG
        
        kf = karamul(k, lf)
        kg = karamul(k, lg)
        
        kf = [elt << (sc1 - sc) for elt in kf]
        kg = [elt << (sc1 - sc) for elt in kg]
        
        lF = [lF[i] - kf[i] for i in range(len(kf))]
        lG = [lG[i] - kg[i] for i in range(len(kg))]
    return lF, lG

def FG_step2(lf, lg, n):
    l = 1
    for k in range(logn[n], 0, -1):
        temp_l = lift(lf[k], l)
        temp_c = conjugate (lg[k - 1], 2 * l)
        lF = karamul(temp_l, temp_c)
        
        temp_l = lift(lg[k], l)
        temp_c = conjugate(lf[k - 1], 2 * l)
        lG = karamul(temp_l, temp_c)
        
        lF,lG = reduce(lF,lG,lf[k - 1],lg [k - 1],2 * l)
        l *= 2
        
        lf [k - 1]= [lF[i] for i in range(l)]
        lg [k - 1]= [lG[i] for i in range(l)]
    return lF, lG

def check_keys(f, g, F, G):
    success=True
    fG = karamul(f, G)
    gF = karamul(g, F)
    
    fGgF=[0 for i in range(len(fG))]
    for i in range(len(fG)):
        fGgF[i] = fG[i] - gF[i]
    if fGgF [0]!=q:
        success=False
    for i in range(1, len(fGgF)):
        if fGgF[i]!=0:
            success=False
    return success

def gen_FG(f,g, n):
    f_copy = deepcopy(f)
    g_copy = deepcopy(g)
    
    lf, lg = FG_step1(f, g, len(f))
    
    gcd,ls,lt = xgcd(lf[logn[n]][0], lg[logn[n]][0])
    success = True
    if gcd!=1:
        return False,[0],[0]
    
    lf[logn[n]][0] = -lt * q
    lg[logn[n]][0] = ls * q
    
    lF, lG = FG_step2 (lf, lg, n)
    
    success = check_keys(f_copy, g_copy, lF, lG)
    return success, lF, lG

def keygen(n):
    success = False
    i = 0
    while success==False:
        f, g, h = gen_fgh(n)
        success, F, G = gen_FG(f, g, n)  
    return f, g, h, F, G