#############################################################################
##
#W  kronecker_graeffe.g           Andrés Herrera-Poyatos <andreshp9@gmail.com>
#W                                Pedro A. Garcia-Sanchez <pedro@ugr.es>
##
## Algorithm to check whether a polynomial is Kronecker proposed in [1].
##
#############################################################################

########################################################################
##
#F IsKroneckerPolynomialGraeffe(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials.
##
##   The algorithm is an improved version of an algorithm proposed in [1]. 
##   It uses the Graeffe method, which, given a polynomial f, provides 
##   another polynomial whose roots are the squares of the roots of f.
##
#########################################################################
IsKroneckerPolynomialGraeffe := function(f)
    local x, sf, is_kronecker;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    if IsZero(f) then
        return false;
    fi;
    
    x:=IndeterminateOfLaurentPolynomial(f);

    # Take the largest square free divisor.
    sf:=f/Gcd(f,Derivative(f));

    # Remove the factor x.
    if Value(sf, 0) = 0 then
        sf := sf / x;
    fi;
    
    # Make sf a primitive polynomial
    sf := sf / Gcd(CoefficientsOfUnivariatePolynomial(sf));
    
    # Check whether sf is monic or not
    if LeadingCoefficient(sf) <> 1 then
        return false;
    fi;    
        
    is_kronecker := function(f)
        local f1, qs, fs, fn, fp, factors_graeffe;
        
        # Remove the factors x-1 and x+1
        if Value(f, 1) = 0 then
            f := f/(x-1);
        fi;
        if Value(f, -1) = 0 then
            f := f/(x+1);
        fi;
        # Check if the polynomial is a constant.
        if Degree(f) = 0 then
            return true;
        fi;
    
        # Check if the polynomial has even degree and is self reciprocal.
        if Degree(f) mod 2 <> 0 or not(IsSelfReciprocalUnivariatePolynomial(f)) then
            return false;
        fi;
        
        # Apply the graeffe method to f. Note that if f(+-x) = f1, 
        # then f is Kronecker.
        f1 := GraeffePolynomial(f);
        if f = f1 then
            return true;
        elif Value(f,-x) = f1 then
            return true;
        fi;
        
        # Assume that f is Kronecker.
        # Find the decomposition of f = fs*fp*fn, where:
        # - fs is the part of f which verifies Graeffe(fs) = (qs(x))^2
        # - fp is the part of f which verifies Graeffe(fp) = fp.
        # - fn is the part of f which verifies Graeffe(fn)(x) = fn(-x).
        qs := Gcd(f1,Derivative(f1));        
        fs := Value(qs, x^2);
        fp := Gcd(f / fs, f1);
        fn := f / (fs * fp);
                
        if fp = f or fn = f then
            # We must have f = fp or f = fn, but we obtained Graeffe(f) != f, 
            # Graeffe(f) != f(-x), a contradiction. Hence f is not Kronecker.
            return false;
        else
            # f is Kronecker if and only if Graeffe(fp) = fp, Graeffe(fn) = (-1)^{deg fn}fn(-x)
            # and qs is Kronecker.
            if GraeffePolynomial(fp) <> fp then
                return false;
            fi;
            if fn <> x+1 and GraeffePolynomial(fn) <> Value(fn, -x) then
                return false;
            fi;            
            if GraeffePolynomial(fs) <> qs^2 then
                return false;
            fi;
            
            return is_kronecker(qs);
        fi;
    end;
    
    # Apply the Graeffe algorithm for square free polynomials
    return is_kronecker(sf);    
end;

########################################################################
##
#F IsKroneckerPolynomialGraeffeOriginal(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials.
##
##   The algorithm was proposed in [1]. It uses the Graeffe method,
##   which, given a polynomial f, provides another polynomial whose roots
##   are the squares of the roots of f.
##
#########################################################################
IsKroneckerPolynomialGraeffeOriginal := function(f)
    local x, sf, is_kronecker;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    if IsZero(f) then
        return false;
    fi;
    
    x:=IndeterminateOfLaurentPolynomial(f);

    # Take the largest square free divisor.
    sf:=f/Gcd(f,Derivative(f));

    # Remove the factor x.
    if Value(sf, 0) = 0 then
        sf := sf / x;
    fi;
    
    # Make sf a primitive polynomial
    sf := sf / Gcd(CoefficientsOfUnivariatePolynomial(sf));
    
    # Check whether sf is monic or not
    if LeadingCoefficient(sf) <> 1 then
        return false;
    fi;    
        
    is_kronecker := function(f)
        local f1, qs, fs, fn, fp, factors_graeffe;
        
        # Remove the factors x-1 and x+1
        if Value(f, 1) = 0 then
            f := f/(x-1);
        fi;
        if Value(f, -1) = 0 then
            f := f/(x+1);
        fi;
        # Check if the polynomial is a constant.
        if Degree(f) = 0 then
            return true;
        fi;
    
        # Check if the polynomial has even degree and is self reciprocal.
        if Degree(f) mod 2 <> 0 or not(IsSelfReciprocalUnivariatePolynomial(f)) then
            return false;
        fi;
        
        # Apply the graeffe method to f. Note that if f(+-x) = f1, 
        # then f is Kronecker.
        f1 := GraeffePolynomial(f);
        if f = f1 then
            return true;
        elif Value(f,-x) = f1 then
            return true;
        fi;
        
        # Assume that f is Kronecker.
        # Find the decomposition of f = fs*fp*fn, where:
        # - fs is the part of f which verifies Graeffe(fs) = (qs(x))^2
        # - fp is the part of f which verifies Graeffe(fp) = fp.
        # - fn is the part of f which verifies Graeffe(fn)(x) = fn(-x).
        qs := Gcd(f1,Derivative(f1));        
        fs := Value(qs, x^2);
        fp := Gcd(f / fs, f1);
        fn := f / (fs * fp);
                
        if fp = f or fn = f then
            # We must have f = fp or f = fn, but we obtained Graeffe(f) != f, 
            # Graeffe(f) != f(-x), a contradiction. Hence f is not Kronecker.
            return false;
        else
            # f is Kronecker if and only if Graeffe(fp) = fp, Graeffe(fn) = (-1)^{deg fn}fn(-x)
            # and qs is Kronecker.
            if GraeffePolynomial(fp) <> fp then
                return false;
            fi;
            if fn <> x+1 and GraeffePolynomial(fn) <> Value(fn, -x) then
                return false;
            fi;            
            if GraeffePolynomial(fs) <> qs^2 then
                return false;
            fi;
            
            return is_kronecker(qs);
        fi;
    end;
    
    # Apply the Graeffe algorithm for square free polynomials
    return is_kronecker(sf);    
end;
