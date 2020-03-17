import math

def decimalToVector(i, l):
    vec = [0] * l
    c = 0

    while i > 0 and c < len(vec):
        vec[-(c + 1)] = i % 2
        i //= 2
        c += 1

    return vec

def vectorToDecimal(i):
    c = 0
    p = 0
    for j in reversed(i):
        c += j * 2**p
        p += 1
    return c

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

def calculate_r(rho, func):
    r = 1
    kr = 0
    flag = True
    while(flag or kr < rho):
        r += 1
        kr = func(r)
        if(kr == rho):
            return r
        flag = False
    return False

def calculate_r_encode(rho):
    return calculate_r(rho, lambda r: 2**r - r - 1)

def calculate_r_decode(rho):
    return calculate_r(rho, lambda r: 2**r - 1)

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

def list_dot(a, b):
    c = 0
    for i in range(len(a)):
        c += a[i] * b[i]
    return c

def is_zero(a):
    for i in a:
        if i != 0:
            return False
    return True

def hammingDecoder(v):
    n = len(v)
    v = v.copy()
    r = calculate_r_decode(n)

    if(not(r)):
        return []
    
    H = []

    # generates H of correct size
    for i in range(r):
        H.append(['-'] * (2**r - 1))

    #fills H
    for i in range(1, 2**r):
        vec = decimalToVector(i, r)
        for j in range(len(H)):
            H[j][i - 1] = vec[j]

    dot = []
    
    #dots rows of H with v
    for i in H:
        dot.append(list_dot(i, v) % 2)

    if(is_zero(dot)):
        return v
    
    bi = vectorToDecimal(dot) - 1

    if(v[bi] == 0):
        v[bi] = 1
    else:
        v[bi] = 0
    
    return v

def messageFromCodeword(c):
    n = len(c)
    r = calculate_r_decode(n)
    if(not(r)):
        return []
    
    m = []

    for i in range(n):
        if(not(math.log(i + 1, 2).is_integer())):
            m.append(c[i])
            
    return m

def message(a):
    rho = len(a)
    r = 1
    flag = True
    while(flag):
        r += 1
        k = 2**r - r - 1
        kr = (2**r - 2 * r - 1)
        flag = not((kr  == (k - r)) and kr >= rho)
    # m = [rho in r bits][data][padding zeroes]
    m = decimalToVector(rho, r) + a + [0]*(k - (r + rho))
    return m

def dataFromMessage(m):
    n = len(m)
    r = calculate_r_encode(n)

    if(not(r)):
        return []
    
    rho = vectorToDecimal(m[0:r])
    
    if(n - r < rho):
        return []

    return m[r:r + rho]

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