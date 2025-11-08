# EDS-DLP example over a small finite field
# Curve: y^2 = x^3 + A*x + B mod p

p = 23        # Finite field
A = 1
B = 1

# Base point P and Q = k*P
P = (3, 10)
k = 2  # secret scalar
Q = (7, 12)  # precomputed for small example

# Modular inverse
def modinv(a, p):
    return pow(a, -1, p)

# Point addition 
def point_add(P, Q):
    if P == Q:
        x, y = P
        lam = ((3*x*x + A) * modinv(2*y, p)) % p
    else:
        x1, y1 = P
        x2, y2 = Q
        lam = ((y2 - y1) * modinv(x2 - x1, p)) % p
    x3 = (lam*lam - P[0] - Q[0]) % p
    y3 = (lam*(P[0] - x3) - P[1]) % p
    return (x3, y3)

# Scalar multiplication
def scalar_mult(n, P):
    R = P
    for _ in range(n-1):
        R = point_add(R, P)
    return R

# Recursive Division polynomials
psi_cache = {}

def psi(n, P):
    x, y = P
    key = (n, x, y)
    if key in psi_cache:
        return psi_cache[key]

    if n == 0:
        val = 0
    elif n == 1:
        val = 1
    elif n == 2:
        val = (2*y) % p
    elif n == 3:
        val = (3*pow(x,4,p) + 6*A*pow(x,2,p) + 12*B*x - A*A) % p
    elif n == 4:
        val = (4*y*(pow(x,6,p) + 5*A*pow(x,4,p) + 20*B*pow(x,3,p) - 5*A*A*pow(x,2,p) - 4*A*B*x - 8*B*B - A**3)) % p
    else:
        # Recurrence: W_{m+n}*W_{m-n} = W_{m+1}*W_{m-1}*W_n^2 - W_{n+1}*W_{n-1}*W_m^2
        if n % 2 == 1:
            k = (n-1)//2
            val = (psi(k+2, P)*pow(psi(k, P),3,p) - psi(k-1,P)*pow(psi(k+1,P),3,p)) % p
        else:
            k = n//2
            psi2 = psi(2, P)
            numerator = (psi(k+2,P)*pow(psi(k-1,P),2,p) - psi(k-2,P)*pow(psi(k+1,P),2,p)) % p
            val = (psi(k,P) * numerator * modinv(psi2,p)) % p

    psi_cache[key] = val
    return val

# Compute W_n(P) and W_n(Q)
n_max = 10
W_P = [psi(n, P) for n in range(1, n_max+1)]
W_Q = [psi(n, Q) for n in range(1, n_max+1)]

print("n | W_n(P) | W_n(Q)")
for i in range(n_max):
    print(f"{i+1} | {W_P[i]}      | {W_Q[i]}")

# Check EDS-DLP identity: W_n(Q) * W_k(P)^{n^2} == W_{n*k}(P) mod p
print("\nChecking EDS-DLP identity with k =", k)
for n in range(1, n_max+1):
    WnkP = psi(n*k, P)
    WkP = psi(k, P)
    lhs = (W_Q[n-1] * pow(WkP, n**2, p)) % p
    print(f"n={n}: LHS={lhs} | W_{n*k}(P)={WnkP} | Match: {lhs==WnkP}")

