#############################################################################
##
#W  kronecker_pieter.g           Andrés Herrera-Poyatos <andreshp9@gmail.com>
#W                              Pedro A. Garcia-Sanchez <pedro@ugr.es>
##
## Algorithm to check whether a polynomial is Kronecker proposed by Pieter Moore.
##
#############################################################################

########################################################################
##
#F IsKroneckerPolynomialPieter(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials.
##
#########################################################################
#DeclareGlobalFunction("IsKroneckerPieter");

########################################################################
##
#F IsKroneckerPolynomialPieter(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials.
##
##   The algorithm was proposed by Pieter Moree. It uses an upper bound onk
##   Degree(f) in terms of f''(1) / f(1) for Kronecker polynomials.
##   The bound was established by A. Herrera-Poyatos and P. Moree.
##   First, we compute the greatest square free factor of f, denoted sf.
##   Now iteratively we compute sf:= sf/gcd(sf, x^n -1). If sf is Kronecker,
##   then at some point we obtain Degree(sf) = 0. Otherwise, Bound(sf)
##   is surpassed by n at some point.
##
#########################################################################
IsKroneckerPolynomialPieter := function(f)
    local sf, d, gd, x, Bound;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    if IsZero(f) then
        return false;
    fi;
        
    x:=IndeterminateOfLaurentPolynomial(f);
    
    # Take the largest square free divisor.
    sf := f / Gcd(f, Derivative(f));

    # Remove the factors x, x+1 and x-1.
    if Value(sf, 0) = 0 then
        sf := sf / x;
    fi;
    if Value(sf, 1) = 0 then
        sf := sf/(x-1);
    fi;
    if Value(sf, -1) = 0 then
        sf := sf/(x+1);
    fi;
    
    # Check if the polynomial is constant.
    if Degree(sf) = 0 then
        return true;
    fi;

    # Check if the polynomial has even degree and is self reciprocal.
    if Degree(sf) mod 2 <> 0 or not(IsSelfReciprocalUnivariatePolynomial(sf)) then
        return false;
    fi;
    
    #########################################
    # Bound(p)
    # Bound for Pieter's algorithm.
    # For Kronecker polynomials it is verified
    # Bound(p) >= min{\Psi(j): \Phi_j divides p}.
    #########################################
    Bound := function(p)
        local d, p21, p1;

        d := DegreeOfUnivariateLaurentPolynomial(p);
        p21 := Value(Derivative(Derivative(p)), 1);
        p1 := Value(p,1);
        return CeilingOfRational((12*p21/p1-3*d^2+6*d)/d);
    end;
    
    d:=3;
    while true do
        if Bound(sf) < d+1 then
            return false;
        fi;
        gd := Gcd(sf, x^d-1);
        if sf = gd then
            return true;
        fi;
        d:=d+1;
        sf:=sf/gd;
    od;
end;

