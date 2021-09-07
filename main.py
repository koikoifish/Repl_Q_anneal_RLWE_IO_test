import numpy as np
import os

def crange(coeffs, q):
    coeffs = np.where((coeffs >= 0) & (coeffs <= q // 2), coeffs, coeffs - q)

    return coeffs


class Rq(object):
    '''
    Ring-Polynomial: Fq[x] / (x^n + 1)
        range of the reminder is set to (−q/2, q/2]
    '''
    def __init__(self, coeffs, q):
        '''
        # Args
            coeffs: coefficients array of a polynomial
            q: modulus
        '''
        n = len(coeffs)  # degree of a polynomial

        f = np.zeros((n + 1), dtype=np.int64)  # x^n + 1
        f[0] = f[-1] = 1
        f = np.poly1d(f)
        self.f = f

        self.q = q
        coeffs = np.array(coeffs, dtype=np.int64) % q
        coeffs = crange(coeffs, q)
        self.poly = np.poly1d(np.array(coeffs, dtype=np.int64))

    def __repr__(self):
        template = 'Rq: {} (mod {}), reminder range: ({}, {}]'
        return template.format(self.poly.__repr__(), self.q, -self.q // 2,
                               self.q // 2)

    def __len__(self):
        return len(self.poly)  # degree of a polynomial

    def __add__(self, other):
        coeffs = np.polyadd(self.poly, other.poly).coeffs
        return Rq(coeffs, self.q)

    def __mul__(self, other):
        q, r = np.polydiv(np.polymul(self.poly, other.poly), self.f)
        coeffs = r.coeffs
        return Rq(coeffs, self.q)

    def __rmul__(self, integer):
        coeffs = (self.poly.coeffs * integer)
        return Rq(coeffs, self.q)

    def __pow__(self, integer):
        if integer == 0:
            return Rq([1], self.q)
        ret = self
        for i in range(integer - 1):
            ret *= ret
        return ret


class RLWE:
    def __init__(self, n, p, t, std):
        assert np.log2(n) == int(np.log2(n))
        self.n = n
        self.p = p
        self.t = t
        self.std = std

    def generate_keys(self):
        s = discrete_gaussian(self.n, self.p, std=self.std)
        e = discrete_gaussian(self.n, self.p, std=self.std)

        a1 = discrete_uniform(self.n, self.p)
        a0 = -1 * (a1 * s + e)
        return (s, (a0, a1))  # (secret, public)

    def encrypt(self, m, a):
        '''
        # Args:
            m: plaintext (mod t)
            a: public key (a0, a1)
        '''
        a0, a1 = a
        e = [discrete_gaussian(self.n, self.p, std=self.std) for _ in range(3)]

        m = Rq(m.poly.coeffs, self.p)

        return (a1 * e[0], np.ceil(self.t / 2) * m + a0 * e[0])

    def decrypt(self, c, s):
        '''
        # Args:
            c: ciphertext (c0, c1, ..., ck)
            s: secret key
        '''

        cc0, cc1 = c
        u = cc0
        v = cc1
        u = Rq(u.poly.coeffs, self.t)
        v = Rq(v.poly.coeffs, self.t)
        s = Rq(s.poly.coeffs, self.t)
        m = v + u * s
        m = Rq(m.poly.coeffs, self.t)
        f = m.poly.coeffs
        out = Rq(np.ones(n), t)
        coeffs = out.poly.coeffs

        for i in range(len(f)):
            if (f[i] >= np.ceil(self.t / 4) or f[i] <= -np.ceil(self.t / 4)):
                coeffs[i] = 1
            else:
                coeffs[i] = 0
        out = Rq(out.poly.coeffs, self.t)
        return (out)

    def add(self, c0, c1):
        '''
        # Args:
            c0: ciphertext (c0, c1, ..., ck)
            c1: ciphertext (c'0, c'1, ..., c'k')
        '''
        c = ()

        k0 = len(c0)  # not necessary to compute (len - 1)
        k1 = len(c1)

        if k0 > k1:
            (c0, c1) = (c1, c0)  # c0 is always shorter

        for _ in range(abs(k0 - k1)):
            c0 += (Rq([0], self.p), )  # add 0 to shorter ciphertext

        for i in range(len(c0)):
            c += (c0[i] + c1[i], )

        return c

    def mul(self, c0, c1):
        '''
        # Args:
            c0: ciphertext (c0, c1, ..., ck)
            c1: ciphertext (c'0, c'1, ..., c'k')
        '''
        c = ()

        k0 = len(c0) - 1
        k1 = len(c1) - 1

        for _ in range(k1):
            c0 += (Rq([0], self.p), )

        for _ in range(k0):
            c1 += (Rq([0], self.p), )

        for i in range(k0 + k1 + 1):
            _c = Rq([0], self.p)
            for j in range(i + 1):
                _c += c0[j] * c1[i - j]
            c += (_c, )

        return c


def discrete_gaussian(n, q, mean=0., std=1.):
    coeffs = np.round(std * np.random.randn(n))
    return Rq(coeffs, q)


def discrete_uniform(n, q, min=0., max=None):
    if max is None:
        max = q
    coeffs = np.random.randint(min, max, size=n)
    return Rq(coeffs, q)


if __name__ == '__main__':
    n = 8  # power of 2
    q = 49  # prime number, q = 1 (mod 2n)
    t = q  
    std = 1  # standard deviation of Gaussian distribution
    #fail =0
    rlwe = RLWE(n, q, t, std)
    #for i in range(10):
      #rlwe = RLWE(n, q, t, std)
      #(sec, pub) = rlwe.generate_keys()
      #print(sec)
      #print(pub)

      #m0 = Rq(np.random.randint(2, size=n), t)

      #c0 = rlwe.encrypt(m0, pub)

      #m_0 = rlwe.decrypt(c0, sec)

      #print(m0)
      #print(c0)
      #print(m_0)
      #if (list(m_0.poly.coeffs) != list(m0.poly.coeffs)):
          #print(1)
          #fail =1
          #print(m0)
          #print(m_0)
    #print(fail)    
    #print()

    df=np.genfromtxt("Shortvects.txt",dtype=int,delimiter=',',usecols=(4,5,6,7,8,9,10,11))
    bigarr = df
    for i in range(len(bigarr)):
      arr = bigarr[i]
      #write something to eat in the csv file?
      testkey = Rq(arr, t)
      pub0 = [22,   4,  -9,  -1,  24, -16,   3,  -7] #this is a 
      pub0 = Rq(pub0,t)
      pub1 = [22, -16, 19, -12, -20, -12, -16, 18] #this is -b
      pub1 = Rq(pub1,t)

      qual=(pub0*testkey + pub1).poly.coeffs #as-b ~ -e
      badentry = 0
      for j in range(len(qual)):
          if(qual[j] >= np.ceil(q/4)) or (qual[j] <= -np.ceil(q/4)):
            badentry += 1
      if(badentry <=1):
        #print(testkey)
        #print(qual)
        pub = (pub1,pub0)
        m0 = Rq(np.random.randint(2, size=n), t)

        c0 = rlwe.encrypt(m0, pub)
        testkey = Rq([-1, -1,  0, -1,  0,  0, -1, -1],t)
        m_0 = rlwe.decrypt(c0, testkey)
        if (list(m_0.poly.coeffs) == list(m0.poly.coeffs)):
          print(i)
          print(m0)
          print(m_0)
