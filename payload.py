import math

def decimalToVector(i, l):
    vec = [0] * l
    c = 0

    while i > 0 and c < len(vec):
        vec[-(c + 1)] = i % 2
        i //= 2
        c += 1

    return vec

def repetitionEncoder(m, n):
    return m * n

def repetitionDecoder(v):
    c = 0
    n = len(v)
    for i in v:
        if(i == 1):
            c += 1
    
    if(c > n // 2):
        return [1]
    elif(n - c > n // 2):
        return [0]
    
    return []

def calculate_r_encode(rho):
    r = 1
    kr = 0
    flag = True
    while(flag or kr < rho):
        r += 1
        kr = (2**r - r - 1)
        if(kr == rho):
            return r
        flag = False
    return False

def calculate_r_decode(rho):
    r = 1
    kr = 0
    flag = True
    while(flag or kr < rho):
        r += 1
        kr = (2**r - 1)
        if(kr == rho):
            return r
        flag = False
    return False

def message(a):
    rho = len(a)
    r = 1
    flag = True
    while(flag):
        r += 1
        k = 2**r - r - 1
        kr = (2**r - 2 * r - 1)
        flag = not((kr  == (k - r)) and kr >= rho)
    m = decimalToVector(rho, r) + a + [0]*(k - (r + rho))
    return m

def list_and(a, b):
    m = []
    for i in range(len(a)):
        m.append(a[i] & b[i])
    return m

def count_ones(a):
    c = 0
    for i in a:
        if(i == 1):
            c += 1
    return c

def hammingEncoder(m):
    k = len(m)
    r = calculate_r_encode(k)  
    if(not(r)):
        return []

    G = hammingGeneratorMatrix(r)
    Gt = list(zip(*G))
    l = len(G[0])
    out = [0] * l
    
    for i in range(len(Gt)):
        mask = Gt[i]
        co = count_ones(list_and(m, mask))
        if(math.log(i + 1, 2).is_integer()):
            if(co % 2 == 0):
                out[i] = 0
            else:
                out[i] = 1
        else:
            out[i] = co
        
    return out

def hammingDecoder(v):
    n = len(v)
    r = calculate_r_decode(n)

    if(not(r)):
        return []
    
    Ht = []
    for i in reversed(range(1, 2**r)):
        Ht.append(decimalToVector(i, 3))
    
    for i in Ht:
        print(i)

    return []

def messageFromCodeword(c):
    return []

def dataFromMessage(m):
    return []

#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1, 2^r-r-1, 3)-Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G

hammingDecoder([1, 1, 1, 0, 0, 0, 0])