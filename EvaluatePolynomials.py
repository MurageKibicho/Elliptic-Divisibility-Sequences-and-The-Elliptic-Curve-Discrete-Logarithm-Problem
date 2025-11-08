class EllipticCurve:
    def __init__(self, a1, a2, a3, a4, a6, p):
        # For curve: y^2 + a1*xy + a3*y = x^3 + a2*x^2 + a4*x + a6
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.a4 = a4
        self.a6 = a6
        self.p = p

class DivisionPolynomials:
    def __init__(self, curve):
        self.curve = curve
        self.p = curve.p
        
    def psi(self, n, x, y):
        """Compute division polynomial ψ_n at point (x,y) for the specific curve"""
        p = self.p
        a1, a2, a3, a4, a6 = self.curve.a1, self.curve.a2, self.curve.a3, self.curve.a4, self.curve.a6
        
        if n == 0:
            return 0
        elif n == 1:
            return 1
        elif n == 2:
            return (2*y + a1*x + a3) % p
        elif n == 3:
            # ψ_3 = 3x^4 + b2*x^3 + 3b4*x^2 + 3b6*x + b8
            b2 = (a1*a1 + 4*a2) % p
            b4 = (a1*a3 + 2*a4) % p
            b6 = (a3*a3 + 4*a6) % p
            b8 = (a1*a1*a6 - a1*a3*a4 + 4*a2*a6 + a2*a3*a3 - a4*a4) % p
            
            term1 = (3 * pow(x, 4, p)) % p
            term2 = (b2 * pow(x, 3, p)) % p
            term3 = (3 * b4 * x*x) % p
            term4 = (3 * b6 * x) % p
            term5 = b8
            
            return (term1 + term2 + term3 + term4 + term5) % p
        elif n == 4:
            # ψ_4 = (2x^6 + b2*x^5 + 5b4*x^4 + 10b6*x^3 + 10b8*x^2 + (b2*b8 - b4*b6)x + b4*b8 - b6^2) * ψ_2
            b2 = (a1*a1 + 4*a2) % p
            b4 = (a1*a3 + 2*a4) % p
            b6 = (a3*a3 + 4*a6) % p
            b8 = (a1*a1*a6 - a1*a3*a4 + 4*a2*a6 + a2*a3*a3 - a4*a4) % p
            
            poly = (2 * pow(x, 6, p)) % p
            poly = (poly + b2 * pow(x, 5, p)) % p
            poly = (poly + 5 * b4 * pow(x, 4, p)) % p
            poly = (poly + 10 * b6 * pow(x, 3, p)) % p
            poly = (poly + 10 * b8 * x*x) % p
            poly = (poly + (b2*b8 - b4*b6) * x) % p
            poly = (poly + (b4*b8 - b6*b6)) % p
            
            psi_2 = self.psi(2, x, y)
            return (poly * psi_2) % p
        else:
            # Use recurrence for larger n
            if n % 2 == 1:  # odd
                k = (n - 1) // 2
                psi_k = self.psi(k, x, y)
                psi_k1 = self.psi(k+1, x, y)
                psi_k2 = self.psi(k+2, x, y)
                psi_km1 = self.psi(k-1, x, y)
                
                term1 = (psi_k2 * pow(psi_k, 3, p)) % p
                term2 = (psi_km1 * pow(psi_k1, 3, p)) % p
                return (term1 - term2) % p
            else:  # even
                k = n // 2
                psi_k = self.psi(k, x, y)
                psi_k1 = self.psi(k+1, x, y)
                psi_k2 = self.psi(k+2, x, y)
                psi_km1 = self.psi(k-1, x, y)
                psi_km2 = self.psi(k-2, x, y)
                psi_2 = self.psi(2, x, y)
                
                numerator = (psi_k2 * pow(psi_km1, 2, p) - psi_km2 * pow(psi_k1, 2, p)) % p
                denominator_inv = pow(psi_2, p-2, p)
                return (psi_k * numerator * denominator_inv) % p

def main():
    # Curve from paper: y^2 + xy + y = x^3 + x^2 + 21x over F_23
    # Convert to: y^2 + 1*xy + 1*y = x^3 + 1*x^2 + 21*x + 0
    p = 23
    curve = EllipticCurve(a1=1, a2=1, a3=1, a4=21, a6=0, p=p)
    div_poly = DivisionPolynomials(curve)
    
    print("=== Computing Division Polynomial Sequences ===")
    print("Curve: y^2 + xy + y = x^3 + x^2 + 21x over F_23")
    print()
    
    # Points from the paper
    P = (0, 0)      # Point P
    Q = (18, 14)    # Point Q = [k]P  
    QP = (21, 0)    # Point Q + P
    
    print("Computing sequence for ψₙ(P) at P = (0,0):")
    psi_P_sequence = []
    for n in range(25):  # Compute first 10 terms
        val = div_poly.psi(n, P[0], P[1])
        psi_P_sequence.append(val)
        print(f"ψ_{n} = {val}")
    
    print("\nComputing sequence for ψₙ(Q) at Q = (18,14):")
    psi_Q_sequence = []
    for n in range(25):
        val = div_poly.psi(n, Q[0], Q[1])
        psi_Q_sequence.append(val)
        print(f"ψ_{n} = {val}")
    
    print("\nComputing sequence for ψₙ(Q+P) at Q+P = (21,0):")
    psi_QP_sequence = []
    for n in range(25):
        val = div_poly.psi(n, QP[0], QP[1])
        psi_QP_sequence.append(val)
        print(f"ψ_{n} = {val}")
    
    # Verify against paper values
    print("\n=== Verification against paper ===")
    print("Paper gives:")
    print("ψₙ(P): 0, 1, 1, 22, 2, ...")
    print(f"Our computed: {psi_P_sequence}")
    
    print("ψₙ(Q): 0, 1, 1, 20, 1, ...")
    print(f"Our computed: {psi_Q_sequence}")
    
    print("ψₙ(Q+P): 0, 1, 22, 11, 18, ...")
    print(f"Our computed: {psi_QP_sequence}")

if __name__ == "__main__":
    main()
