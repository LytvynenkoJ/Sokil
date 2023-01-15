from sokil_constants import logn, FFT_TABL

def cpoly_fft(f,n):
    log2n = logn[n]
    hn = n // 2
    t = hn
    u = 1
    m = 2
    while u < log2n:
        ht = t // 2; hm = m // 2
        i1 = 0; j1 = 0
        while i1 < hm:
            j2 = j1 + ht
            s_re = FFT_TABL [(m + i1) * 2]
            s_im = FFT_TABL [(m + i1) * 2 + 1]
            j = j1
            while j<j2:
                x_re = f[j]; x_im = f[j + hn]
                y_re = f[j + ht]; y_im = f[j + ht + hn]
                y_re, y_im = s_re * y_re - s_im * y_im, s_re * y_im + s_im * y_re     # CMUL(y_re, y_im, s_re, s_im)
                f[j], f[j + hn] = x_re + y_re, x_im + y_im                            # CADD (x_re, x_im, y_re, y_im)
                f[j + ht], f[j + ht + hn] = x_re - y_re, x_im - y_im                  # CSUB( x_re, x_im, y_re, y_im)
                j = j + 1
            i1 = i1 + 1
            j1 = j1 + t
        t = ht
        u = u + 1
        m = m * 2
    return f

def cpoly_inv_fft(f,n):
    log2n = logn[n]; t = 1; m = n; hn = n // 2; u = log2n
    while u > 1:
        hm = m // 2
        dt = t * 2
        i1 = 0; j1 = 0
        while j1 < hn:
            j2 = j1 + t
            s_re = FFT_TABL[(hm + i1) * 2]
            s_im = -FFT_TABL[(hm + i1) * 2 + 1]
            j = j1
            while j < j2:
                x_re = f[j]; x_im = f[j + hn]
                y_re = f[j + t]; y_im = f[j + t + hn]
                f[j], f[j + hn] = x_re + y_re, x_im + y_im                                         # CADD(x_re, x_im, y_re, y_im)
                x_re, x_im = x_re - y_re, x_im - y_im                                              # CSUB( x_re, x_im, y_re, y_im)
                f[j + t], f[j + t + hn] = x_re * s_re - x_im * s_im, x_re * s_im + x_im * s_re     # CMUL0( x_re, x_im, s_re, s_im)
                j = j + 1
            i1= i1 + 1
            j1= j1 + dt
        t = dt; m = hm; u -= 1
    ni=1/(2**(log2n-1))
    for u in range(n):
        f[u]*=ni
    return f

def cpoly_adj(a,n):
    u = n//2
    while u < n:
        a[u] = -a[u]
        u = u + 1
    return a

def сpoly_add_muladj_fft(F,G,f,g, n):
    d=[0 for i in range(n)]
    hn = n // 2
    u = 0
    while u < hn:
        F_re = F[u]; F_im = F[u + hn]
        G_re = G[u]; G_im = G[u + hn]
        f_re = f[u]; f_im = -f[u + hn]
        g_re = g[u]; g_im = -g[u + hn]
        a_re, a_im = F_re * f_re - F_im * f_im, F_re * f_im + F_im * f_re
        b_re, b_im = G_re * g_re - G_im * g_im, G_re * g_im + G_im * g_re
        d[u] = a_re + b_re; d[u + hn] = a_im + b_im
        u = u + 1
    return d

def сpoly_ddiv_fft(a,b):
    n = len(a)
    hn = n // 2
    u = 0
    while u < hn:
        ib = 1 / b [u]
        a[u] *= ib
        a[u + hn] *= ib
        u = u + 1
    return a