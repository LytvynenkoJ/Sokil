from keygen import *
from encoding import encode_i16, decode_i16, encode_i8, decode_i8

class PublicKey:
    """
    This class contains methods for performing public key operations in Sokil.
    """

    def __init__(self, sk, print_format=0):
        """Initialize a public key."""
        self.n = sk.n
        self.h = sk.h
        self.signature_bound = sk.signature_bound
        self.pk_bytelen = sk.pk_bytelen
        self.encoded_pk = ""
        self.encoded_pk_len = 0
        self. print_format = print_format

    def __repr__(self):
        """Print the object in readable form."""
        rep = ""
        if(self.print_format==0 or self.print_format==2):
            rep += "Public key for n = {n}:\n\n".format(n=self.n)
            rep += "h = {h}\n".format(h=self.h)
        if(self.print_format==1 or self.print_format==2):
            if self.encoded_pk!="":
                rep += "\n\nEncoded public key for n = {n}:\n\n".format(n=self.n)
                rep += "Encoded = {enc}\n".format(enc=self.encoded_pk)
            else:
                rep += "Sorry, the key wasn't encoded jet. Use encode() function to encode the key."
        return rep
    
    def setPrintFormat(self, form):
        self.print_format = form%3
        return
    
    def encode(self, h = None):
        if h==None:
            h = self.h
        n = len(h)
        pk = str(hex(0x00 + logn[n]))+" "
        pk_len =1
        l, temp = encode_i16(h)
        pk_len+=l
        pk+=temp
        if h==self.h:
            self.encoded_pk = pk
            self.encoded_pk_len = pk_len
        return pk_len, pk
    
    def decode(self, pk_len = None, pk = None):
        if pk_len==None:
            pk_len = self.encoded_pk_len
        if pk==None:
            pk = self.encoded_pk
        pk2 = pk.split()
        k = str(hex(0x00 + logn[self.n]))
        if pk_len==self.pk_bytelen and pk2[0]==k:
            h = decode_i16(pk2[1:])
            return h
        

class SecretKey:
    """
    This class contains methods for performing
    secret key operations (and also public key operations) in Sokil.
    """

    def __init__(self, n, print_format = 0, polys=None):
        """Initialize a secret key."""
        # Public parameters
        self.n = n
        self.sigma = Params[n]["sigma"]
        self.sigmin = Params[n]["sigmin"]
        self.signature_bound = Params[n]["sig_bound"]
        self.sk_bytelen = Params[n]["sk_bytelen"]
        self.pk_bytelen = Params[n]["pk_bytelen"]
        self.encoded_sk = ""
        self.encoded_sk_len = 0
        self. print_format = print_format

        if polys is None:
            self.f, self.g, self.h, self.F, self.G = keygen(n)
        else:
            [f, g, F, G, h] = polys
            assert all((len(poly) == n) for poly in [f, g, F, G])
            self.f = f[:]
            self.g = g[:]
            self.F = F[:]
            self.G = G[:]
            self.h = h[:]

    def __repr__(self):
        """Print the object in readable form."""
        rep = ""
        if(self.print_format==0 or self.print_format==2):
            rep += "Private key for n = {n}:\n\n".format(n=self.n)
            rep += "f = {f}\n".format(f=self.f)
            rep += "g = {g}\n".format(g=self.g)
            rep += "F = {F}\n".format(F=self.F)
            rep += "G = {G}\n".format(G=self.G)
        if(self.print_format==1 or self.print_format==2):
            if self.encoded_sk!="":
                rep += "\n\nEncoded private key for n = {n}:\n\n".format(n=self.n)
                rep += "Encoded = {enc}\n".format(enc=self.encoded_sk)
            else:
                rep += "Sorry, the key wasn't encoded jet. Use encode() function to encode the key."
        return rep
    
    def setPrintFormat(self, form):
        self.print_format = form%3
        return
    
    def encode(self, f = None, g = None, F = None):
        if f==None:
            f = self.f
        if g == None:
            g = self.g
        if F == None:
            F = self.F
        n = len(f)
        sk = str(hex(0x50 + logn[n]))+" "
        sk_len =1
        l, temp = encode_i8(f, SK_0_BITS)
        sk_len= sk_len + l
        sk = sk + temp
        l, temp = encode_i8(g, SK_0_BITS)
        sk_len= sk_len + l
        sk = sk + temp
        l, temp = encode_i8(F, SK_1_BITS)
        sk_len= sk_len + l
        sk = sk + temp
        if f==self.f and g==self.g and F==self.F:
            self.encoded_sk = sk
            self.encoded_sk_len = sk_len
        return sk_len, sk
    
    def decode(self, sk_len = None, sk = None, n = None):
        if sk_len==None:
            sk_len = self.encoded_sk_len
        if sk==None:
            sk = self.encoded_sk
        sk2 = sk.split()
        k = str(hex(0x50 + logn[self.n]))
        if sk_len==Params[self.n]["sk_bytelen"] and sk2[0]==k:
            l = (self.n*SK_0_BITS)//8+1
            f = decode_i8(sk2[1:l], SK_0_BITS)
            g = decode_i8(sk2[l:2*l-1], SK_0_BITS)
            l = 2*l-1
            F = decode_i8(sk2[l:], SK_1_BITS)
            G = recoveryG(f, g, F)
            return f, g, F, G    