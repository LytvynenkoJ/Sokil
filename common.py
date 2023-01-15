from gmpy2 import mpfr

def split(f):
    """Split a polynomial f in two polynomials.

    Args:
        f: a polynomial

    Format: coefficient
    """
    n = len(f)
    f0 = [f[2 * i + 0] for i in range(n // 2)]
    f1 = [f[2 * i + 1] for i in range(n // 2)]
    return [f0, f1]


def merge(f_list):
    """Merge two polynomials into a single polynomial f.

    Args:
        f_list: a list of polynomials

    Format: coefficient
    """
    f0, f1 = f_list
    n = 2 * len(f0)
    f = [0] * n
    for i in range(n // 2):
        f[2 * i + 0] = f0[i]
        f[2 * i + 1] = f1[i]
    return f


def sqnorm(v):
    """Compute the square euclidean norm of the vector v."""
    res = 0
    for elt in v:
        for coef in elt:
            res += coef ** 2
    return res

def has_Null(f):
    for i in range(len(f)//2):
        if f[i]==0 and f[len(f)//2+i]==0:
            return True
    return False

def poly_mul(x, y):
    t=len(x)
    temp=[mpfr(0) for i in range(2*t)]
    u = 0
    for u in range(t):
        for v in range(t):
            temp[u + v] = temp[u + v] + x[u] * y[v]
    z=[0 for i in range(t)]
    for i in range(t):
        z[i] = temp[i] - temp[i + t]
    return z

def karatsuba(a, b, n):
    if n == 1:
        return [a[0] * b[0], 0]
    else:
        n2 = n // 2
        a0 = a[:n2]
        a1 = a[n2:]
        b0 = b[:n2]
        b1 = b[n2:]
        ax = [a0[i] + a1[i] for i in range(n2)]
        bx = [b0[i] + b1[i] for i in range(n2)]
        a0b0 = karatsuba(a0, b0, n2)
        a1b1 = karatsuba(a1, b1, n2)
        axbx = karatsuba(ax, bx, n2)
        for i in range(n):
            axbx[i] -= (a0b0[i] + a1b1[i])
        ab = [0] * (2 * n)
        for i in range(n):
            ab[i] += a0b0[i]
            ab[i + n] += a1b1[i]
            ab[i + n2] += axbx[i]
        return ab
    
def karamul(a, b):
    n = len(a)
    ab = karatsuba(a, b, n)
    abr = [ab[i] - ab[i + n] for i in range(n)]
    return abr

def ext_gcd(a,b):
    a_ = a; b_ =  b; u = 1; v = 0; s = 0; t = 1; alpha = a; beta = b
    print(a_%2, b_%2)
    if a_%2 == 0 and b_%2 == 0:
        return False, 0, 0
    while a_%2==0:
        a_ = a_ // 2
        if u%2==0 and v%2==0:
            u = u//2; v=v//2
        else:
            u += beta; u = u//2
            v -= alpha; v = v//2
    while a_!=b_:
        if b_%2 == 0:
            b_=b_//2
            if s%2==0 and t%2==0:
                s = s//2; t = t//2
            else:
                s += beta; s = s//2
                t -= alpha; t = t//2
        else:
            if b_<a_:
                b_, a_ = a_, b_
                u, s = s, u
                v, t = t, u
            else:
                b_ -= a_
                s -= u
                t -= v
    gcd = a_
    success = True
    if gcd!=1:
        success = False
    return success, s, t

def xgcd(b, n):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while n != 0:
        q, b, n = b // n, n, b % n
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return b, x0, y0

def lmin(l, t):
    m=0
    for i in range(t):
        if l[i]>m:
            m=l[i]
    temp = m%255
    m = m//255 + int(temp)
    return m

def scale(lf, lg, t):
    f_bits = lmin(lf, t)
    g_bits = lmin(lg, t)
    sc = max(53, f_bits, g_bits) - 53
    return sc

def bitsize(a):
    val = abs(a)
    res = 0
    while val:
        res += 8
        val >>= 8
    return res

def lift(l, t):
    l2=[0 for i in range(2*t)]
    for i in range(t):
        l2[2*i]=l[i]
    return l2

def conjugate(l, t):
    l2 = [0 for i in range(t)]
    for i in range(t):
        if i%2==1:
            l2[i]=-l[i]
        else:
            l2[i] = l[i]
    return l2