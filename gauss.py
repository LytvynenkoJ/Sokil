# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 19:50:14 2023

@author: Yuliia
"""
from kupyna import *
from strumok import *
from sokil_constants import *
import numpy as np

def PsevdoRandomInit(seed):
    temp1=(seed+b' 0x0').decode('utf-8')
    temp2=(seed+b' 0x1').decode('utf-8')
    key = Kupyna(32).hash(temp1)
    iv = Kupyna(32).hash(temp2)
    casted_key = [hex_to_number(number) for number in key]
    casted_iv = [hex_to_number(number) for number in iv]
    strumokInit = Strumok(casted_key, 256, casted_iv)
    return strumokInit

seed=b'0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0 0x0'
strumok=PsevdoRandomInit(seed)

def Rand(strumok_instance, l, buffer):
    if (l-2)%16==0:
        buffer = strumok_instance.update_state()
    z_array=buffer[(l-2)%16:l%16]
    return z_array, buffer

unsigned_part = (1 << 63) - 1

def mkgauss(n):
    g = 1
    if n == 512:
        g= 2
    if n!=2048:
        GAUSS_TABL = tabl_gauss_1024
    else:
        GAUSS_TABL = tabl_gauss_2048
    val = 0
    t = 0
    num=2
    buffer=[0 for i in range(16)]
    while t < g:
        temp, buffer=Rand(strumok, num, buffer)

        r1 = np.int64(temp[0]).item()
        r2 = np.int64(temp[1]).item()
        
        sign = np.sign(r1)
        r1 = r1 & unsigned_part
        r2 = r2 & unsigned_part

        num+=2

        if r1 > GAUSS_TABL[0]:
            k= 1
            while GAUSS_TABL[k] > r2:
                k=k+1
            val = val + sign * k
        t = t + 1
    return int(val)

def poly_mkgauss(n):
    s=0
    f=[0 for i in range(n)]
    for u in range (n-1):
        f[u] = mkgauss(n)
        s = s^f[u]
    success=False
    u=n-1
    while success==False:
        f[u]=mkgauss(n)
        if (s^f[u])%2==1:
            success=True
    return f

def sqnorm(v):
    res = 0
    for elt in v:
        res += elt ** 2
    return res