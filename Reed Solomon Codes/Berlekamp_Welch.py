import numpy as np
import galois
from tqdm import tqdm


# Solve the Ax = b equation using row reduce
def solve(A, b, f):
    A = ff(A)
    b = ff(b).reshape(-1, 1)
    m, n = A.shape

    # Make the Augmented matrix and perform Gaussian Elimination
    Aug = np.hstack([A, b])
    R = ff(Aug).row_reduce()
    lhs = R[:, :n]
    rhs = R[:, n]

    # lhs is the matrix part, rhs is the 'source' (column vector) part

    # Check consistency
    for i in range(m):
        if np.all(lhs[i, :] == 0) and rhs[i] != 0:
            print("No solution")
            exit(1)

    # Identify pivot columns and their rows (non zero values in the RREF of A)
    mp = {}
    for i in range(m):
        nz = np.where(lhs[i, :] != 0)[0]
        if nz.size > 0:
            mp[nz[0]] = i

    pivots = sorted(mp.keys())

    # Construct one solution (free vars set to 0)
    x = ff.Zeros(n)
    for c in pivots:
        x[c] = rhs[mp[c]]

    return x


# Read the file
file_path = 'music4_corrupted.mp3'

with open(file_path, 'rb') as f:
    binary_data = f.read()

binary_data = list(binary_data)
info = []
for i in range(0, len(binary_data), 125):
    a = binary_data[i:i+125] # one alpha
    while len(a) < 125:
        a.append(0)
    info.append(a)

info = np.array(info, dtype=np.int64)

# parameters of the underlying RS code
n, k, m = 125, 100, 8
tc = int(np.floor((n-k)/2))
# finite field
ff = galois.GF(2**m, irreducible_poly="x^8+x^4+x^3+x^2+1")
info = ff(info)

# evaluation set
eval = np.arange(125)
eval = ff(eval)

decoded = []
for x in tqdm(info):
    M = [] # matrix
    z = [] # RHS
    for alpha, y in (zip(eval, x)):
        b = []
        for i in range(tc+k):
            b.append(alpha**i)
        for j in range(tc):
            b.append(-y*(alpha**j))
        M.append(b)
        z.append(y*(alpha**tc))
    # converting to finite field
    M = np.array(M, dtype=np.int64)
    z = np.array(z, dtype=np.int64)
    M = ff(M)
    z = ff(z)

    coeffs = solve(M, z, ff)
    # coefficients of the B and E polynomials
    b_coeffs = coeffs[:tc+k]
    e_coeffs = coeffs[tc+k:]
    e_coeffs = np.concatenate([e_coeffs, [1]])
    B = galois.Poly(b_coeffs, field=ff, order='asc')
    E = galois.Poly(e_coeffs, field=ff, order='asc')
    f, rem = B//E, B%E
    if rem != 0:
        print("B is not divisible by E")
        exit(1)

    for i in range(k):
        decoded.append(f(eval[i]))

decoded = bytes(decoded)
with open("Group4_decoded.mp3", "wb") as f:
    f.write(decoded)