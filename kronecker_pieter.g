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
##   The algorithm was proposed by Pieter Moree. It uses an upper bound on
##   Degree(f) in terms of f''(1) / f(1) for Kronecker polynomials.
##   The bound was established by A. Herrera-Poyatos and P. Moree.
##   First, we compute the greatest square free factor of f, denoted sf.
##   Now iteratively we compute sf:= sf/gcd(sf, x^n -1). If sf is Kronecker,
##   then at some point we obtain Degree(sf) = 0. Otherwise, Bound(sf)
##   is surpassed by n at some point.
##
#########################################################################
IsKroneckerPolynomialPieter := function(f)
    local sf, d, gd, x, Bound, B;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    if IsZero(f) then
        return false;
    fi;
        
    x:=IndeterminateOfLaurentPolynomial(f);
    
    # Check if it is self-reciprocal!!

    
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
    
    # Check if the polynomial is a constant.
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
        return Int(12*p21/p1/d-3*d+6);
    end;
    
    d:=3;
    B:= Bound(sf);
    
    while d <= B do         
        gd := Gcd(sf, x^d-1);
        if Degree(gd) > 0 then
            sf := sf / gd;
            if Degree(sf) = 0 then
                return true;            
            # Check if the polynomial has even degree and is self reciprocal.
            elif Degree(sf) mod 2 <> 0 or not(IsSelfReciprocalUnivariatePolynomial(sf)) then
                return false;
            fi;
            B := Bound(sf);            
        fi;
        d:=d+1;
    od;
    return false;    
end;



########################################################################
##
#F IsKroneckerPolynomialPieterImproved(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials.
##
##   The algorithm was proposed by Pieter Moree. It uses an upper bound on
##   the integers d such that Phi_d divides f.
##   The bound was established by A. Herrera-Poyatos and P. Moree.
##   First, we compute the greatest square free factor of f, denoted sf.
##   Now iteratively we compute sf:= sf/gcd(sf, x^n -1). If sf is Kronecker,
##   then at some point we obtain Degree(sf) = 0. Otherwise, Bound(sf)
##   is surpassed by n at some point.
##
#########################################################################
IsKroneckerPolynomialPieterImproved := function(f)
    local sf, d, gd, x, Bound1, Bound, pi, B;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    if IsZero(f) then
        return false;
    fi;
        
    x:=IndeterminateOfLaurentPolynomial(f);
    
    # Check if it is self-reciprocal!!

    
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
    
    # Check if the polynomial is a constant.
    if Degree(sf) = 0 then
        return true;
    fi;
    
    # Check if the polynomial has even degree and is self reciprocal.
    if Degree(sf) mod 2 <> 0 or not(IsSelfReciprocalUnivariatePolynomial(sf)) then
        return false;
    fi;
    
    #########################################
    # Bound1(p)
    # Bound for Pieter's algorithm.
    # For Kronecker polynomials it is verified
    # Bound1(p) >= min{\Psi(j): \Phi_j divides p}.
    #########################################
    Bound1 := function(p)
        local d, p21, p1;

        d := DegreeOfUnivariateLaurentPolynomial(p);
        p21 := Value(Derivative(Derivative(p)), 1);
        p1 := Value(p,1);
        return 12*p21/p1/d-3*d+6;
    end;
    
    #########################################
    # Bound(p)
    # Bound for the improved Pieter's algorithm.
    # For Kronecker polynomials it is verified
    # Bound(p) >= d
    # for any d such that Phi_d divides p 
    #########################################
    pi := 31416/10000;    
    Bound := function(p)
        local B;
        B := Bound1(p);
        return Minimum(2*Degree(p)^2, Int(B), Int(pi*RootInt(6*B*Degree(p),2)+3));        
    end;
    
    d:=3;
    B:= Bound(sf);
    
    Print(B, sf);
    
    while d <= B do         
        gd := Gcd(sf, x^d-1);
        if Degree(gd) > 0 then
            sf := sf / gd;
            if Degree(sf) = 0 then
                return true;            
            # Check if the polynomial has even degree and is self reciprocal.
            elif Degree(sf) mod 2 <> 0 or not(IsSelfReciprocalUnivariatePolynomial(sf)) then
                return false;
            fi;
            B := Bound(sf);            
        fi;
        d:=d+1;
    od;
    return false;    
end;
