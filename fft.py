# for some computations
import math
import numpy as np

# fast arithmetic (to be improved)
def prod(a, b):

    a = np.remainder(a, S)
    b = np.remainder(b, S)

    return np.remainder(a * b, S)

def pow(b, e):

    return np.remainder(b ** e, S)

def rem(n):

    return np.remainder(n, S)

# whether n is a power of 2
def is_power(n):

    x = int(math.floor(np.log2(n)))
    y = 1 << x
    return (y == n)

# whether or not w is a root of unity of order k
def root_of_uni(w, k):

    # kth power of w is 1
    if not pow(w, k) == 1:
        return False

    # all powers of w < k are not 1
    for i in range(1, k):
        if pow(w, i) == 1:
            return False

    return True;

# split into B (even indeces of A) and C (odd indeces of A)
def split(A):

    B = []
    C = []

    for i in range(0, len(A) / 2):

        B = B + [A[2 * i]]
        C = C + [A[2 * i + 1]]

    return B, C

# hardcoded special cases
def fft_sc(k, A):

    # special case (a)
    if k == 2:
        e0 = A[0] + A[1]
        e1 = A[0] - A[1]
        return [e0, e1]

    # special cases (b), (c), (d)
    # (b) A = (0, ... , 0)      return A
    b = (A[0] == 0)
    # (c) A = (a, 0, ..., 0)    return (a, ... a)
    c = not b
    # (d) A = (a, ...., a)      return (ka, 0, ..., 0)
    d = not b

    for i in range(1, len(A)):

        a = A[i]

        if not a == 0:
            b = False
            c = False

        if not a == A[0]:
            d = False

    if b:
        return A

    if c:
        return [A[0]] * len(A)

    if d:
        end = [0] * (len(A) - 1)
        return [prod(k, A[0])] + end

# k: a power of two
# A: a sequence of length k of elements in S, a finite computation structure
# w: a kth root of unity in S
def fft(k, A, w):

    assert is_power(k)
    assert root_of_uni(w, k)

    # step 0: base cases
    if k == 1:
        return A

    # step .5: hardcoding some special cases for speed up
    # output
    A_ = None
    A_ = fft_sc(k, A)
    if not A_ == None:
        # sol_set[n] = A_
        return A_

    A_ = [None] * len(A)

    # step 1: split
    B, C = split(A)

    # step 2: recurse

    # divide k by 2
    # right shift by n bit is equivalent to division by pow(2, n)
    k = k >> 1
    w_ = pow(w, 2)
    B_ = fft(k, B, w_) #even
    C_ = fft(k, C, w_)

    # step 3:
    for i in range(0, k):
        A_[i] = B_[i] + prod((w ** i), C_[i])
        j = i + k
        A_[j] = B_[i] + prod((w ** j), C_[i])

    return A_

# iterative fft
def it_fft(n, A, w):

    roots = [None] * n
    roots[n - 1] = w

    for i in range(n - 1, 1, -1):
        # roots[i] = roots[i + 1] ** 2
        roots[i - 1] = roots[i] ** 2
        k = 1 << i
        assert(root_of_uni(roots[i], k))

    for s in range(1, n + 1):
        m = 1 << s
        l = m >> 1 # m / 2

        w_ = roots[s - 1]

        for k in range(0, n, m):

            temp = 1

            for j in range(0, l):

                t = prod(w, A[k + j + l])
                u = A[k + j]
                A[k + j] = np.remainder(u + t, S)
                A[k + j + l] = np.remainder(u - t, S)
                temp = prod(w, w_)

    return A


# p: prime of order 256 or 512 bits
# n: a number
# w: 2^nth root of unity
# A: 2^n vector of coeficients fot the polynomial p

p, n, w, A = (5, 2, 2, [1, 4, 0, 0])

# left shift by n bit is equivalent to multiplication by pow(2, n)
k = 1 << n
# finite computation structure: large prime field parameter
S = p

A_ = fft(k, A, w)
print A_
