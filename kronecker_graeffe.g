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
##   The algorithm was proposed in [1]. It uses the Graeffe method,
##   which, given a polynomial f, provides another polynomial whose roots
##   are the squares of the roots of f.
##
#########################################################################
IsKroneckerPolynomialGraeffe := function(f)
    local x, sf, f1, fs, fn, fp, coeffs_fs, factors_graeffe;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    if IsZero(f) then
        return false;
    fi;
    
    x:=IndeterminateOfLaurentPolynomial(f);

    # Take the largest square free divisor.
    sf:=f/Gcd(f,Derivative(f));

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
    
    # Apply the graeffe method to sf. Note that if sf(+-x) = f1, then
    # sf is Kronecker.
    f1:= GraeffePolynomial(sf);
    if sf = f1 then
        return true;
    fi;
    if Value(sf,-x) = f1 then
        return true;
    fi;
    
    # Assume that sf is Kronecker.
    # Find the descomposition of sf = fs*fp*fn, where:
    # - fs is the part of sf which verifies Graeffe(fs) = (g(x))^2
    # - fp is the part of sf which verifies Graeffe(fs) = fs.
    # - fn is the part of sf which vefifies Graeffe(fs)(x) = fs(-x).
    fs := Value(Gcd(Derivative(f1), f1), x^2);
    fp := Gcd(sf / fs, f1);
    fn := sf / (fs * fp);

    factors_graeffe := Difference([fs, fp, fn], [1+0*x]);

    if Length(factors_graeffe) = 1 then
        if IsOne(fs) then
            # We must have f = fp or f = fs, but we obtained Graeffe(sf) != sf, 
            # Graeffe(sf) != sf(-x), a contradiction. Hence sf is not Kronecker.
            return false;
        else
            # sf is Kronecker if and only if f1 / fs(\sqrt{x}) is Kronecker.
            return IsKroneckerPolynomialGraeffe(f1 /  Gcd(f1,Derivative(f1)));
        fi;
    else
        # sf is Kronecker if and only if fs, fp and fn are Kronecker.
        return ForAll(factors_graeffe, IsKroneckerPolynomialGraeffe);
    fi;
end;
