"""SM group theory and fields classes and functions."""

from sympy import Rational

global GUTU1
GUTU1 = True    # Is U(1) coupling normalized for GUT unification?

def CA(N):
    """Quadratic Casimir in adjoint rep C2(A)."""
    if N == 1: 
        return 0  # no self-interaction in U(1)
    else:
        return N
    
def DA(N):
    """Dimension of adjoint rep D(A)."""
    if N == 1: 
        return 1
    else:
        return N**2 - 1


class irrep(object):
    """IRREP of SM group
    
    specified by triplet (D3, D2, Y, multiplicity) of dimensions
    w.r.t. SU(3), SU(2) and hypercharge Y, within
    Q = T3 + Y/2 normalization.
    
    multiplicity -- number of such reps, e.g. n_g
    
    """
    
    def __init__(self, D3, D2, Y, multiplicity=1):
        self.D3 = D3
        self.D2 = D2
        self.Y = Y
        self.multiplicity = multiplicity
        self.D = self.D3 * self.D2  # total dimension
        
    def __str__(self):
        return "{}({},{},{},{})".format(self.__class__.__name__, self.D3, self.D2, self.Y, self.multiplicity)
    
    def C2(self, N):
        """Quadratic Casimir.
        
        By default, for U(1) GUT normalization sqrt(3/5) is included.
        """
        global GUTU1
        if N == 1:
            if GUTU1 == True:
                GUTcsq = Rational(3,5)
            else:
                GUTcsq = 1
            C2 = Rational(1,4)*self.Y**2 * GUTcsq
        elif N == 2:
            C2 = Rational(1,4)*(self.D2**2-1)
        elif N == 3:             
            dim2pq = {1: (0,0), 3: (1,0), 6: (2,0), 8: (1,1)}
            p, q = dim2pq[self.D3]
            C2 = Rational(1,3)*(p**2 + q**3 +3*p + 3*q + p*q)
        return C2 
            
        
    def S2(self, N):
        """Quadratic Dynkin index"""
        return self.multiplicity * self.D * self.C2(N) / DA(N)


def SMDynkin(reps):
    """Total Dynkin index for SM gauge groups, adding all reps."""
    S1, S2, S3 = (0, 0, 0)
    for rep in reps:
        S1 += rep.wS2(1)
        S2 += rep.wS2(2)
        S3 += rep.wS2(3)
    return S1, S2, S3


class Weyl(irrep):
    
    def wS2(self, N):
        """Quadratic Dynkin weighted with n.d.o.f."""
        return Rational(1,2)*self.S2(N)
    

class Dirac(irrep):
    
    def wS2(self, N):
        """Quadratic Dynkin weighted with n.d.o.f."""
        return self.S2(N)

    
class RealScalar(irrep):
    
    def wS2(self, N):
        """Quadratic Dynkin weighted with n.d.o.f."""
        return Rational(1,2)*self.S2(N)
    

class ComplexScalar(irrep):
    
    def wS2(self, N):
        """Quadratic Dynkin weighted with n.d.o.f."""
        return self.S2(N)


