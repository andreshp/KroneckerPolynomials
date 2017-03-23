#############################################################################
##
#W  sturm_kronecker.g           Andrés Herrera-Poyatos <andreshp9@gmail.com>
#W                              Pedro A. Garcia-Sanchez <pedro@ugr.es>
##
## Algorithm to check whether a polynomial is Kronecker based on
## the Sturm sequence.
##
#############################################################################

########################################################################
##
#F SturmSequence(f) returns the Sturm sequence of the polynomial f.
##
#########################################################################
DeclareGlobalFunction("SturmSequence");

########################################################################
##
#F NumberOfRootsOfUnivariatePolynomialInInterval(f, a, b) returns the 
## number of roots that the polynomial f has in the interval [a,b].
## It uses the sturm sequence of the polynomial.
##
#########################################################################
DeclareGlobalFunction("NumberOfRootsOfUnivariatePolynomialInInterval");

########################################################################
##
#F IsKroneckerPolynomialSturm(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials
##
#########################################################################
DeclareGlobalFunction("IsKroneckerPolynomialSturm");

########################################################################
##
#F SturmSequence(f) returns the Sturm sequence of the polynomial f.
##
#########################################################################
SturmSequence := function(p)
    local sst;

    if not(IsUnivariatePolynomial(p)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    sst := [p, Derivative(p)];
    while Degree(sst[Length(sst)]) >= 1 do
        Add(sst, -QuotientRemainder(sst[Length(sst)-1],sst[Length(sst)])[2]);
    od;
    return sst;
end;

########################################################################
##
#F NumberOfRootsOfUnivariatePolynomialInInterval(f, a, b) returns the 
## number of roots that the polynomial f has in the interval [a,b].
## It uses the sturm sequence of the polynomial.
##
#########################################################################
NumberOfRootsOfUnivariatePolynomialInInterval := function(p,a,b)
    local sst, sgnc, valuesa, valuesb, psf;

    if not(IsUnivariatePolynomial(p)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    # Take the largest square free divisor.
    psf:=p/Gcd(p,Derivative(p));

    sgnc := function(numbers)
        local i, n, current_sign, new_sign, numbers_nozeros, count;
        count := 0;
        numbers_nozeros := Filtered(numbers, x -> x <> 0);
        n := Length(numbers_nozeros);
        current_sign := numbers_nozeros[1] > 0;

        for i in [2..n] do
            new_sign := numbers_nozeros[i] > 0;
            if new_sign <> current_sign then
                count := count + 1;
            fi;
            current_sign := new_sign;
        od;
        return count;
    end;

    sst := SturmSequence(p);
    valuesa := List(sst, psf -> Value(psf, a));
    valuesb := List(sst, psf -> Value(psf, b));

    return sgnc(valuesa) - sgnc(valuesb);
end;

########################################################################
##
#F IsKroneckerPolynomialSturm(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials.
##
##   The algorithm was proposed by David Boyd.
##   First we obtain the polynomial sf as the largest square free factor of f.
##   Then we check whether x, x-1 or x+1 divide sf. If it is the case,
##   then we remove these factors from sf. Afterwards if deg sf is 0, then
##   f is Kronecker. If deg sf is not even, then f is not Kronecker.
##   
##   It uses the sturm's algorithm to count roots along with a variable's 
##   change in order to count the roots in the unit circle.
#########################################################################
IsKroneckerPolynomialSturm := function(f)
    if not(IsUnivariatePolynomial(p)) then
        Error("The argument must be a polynomial in one variable");
    fi;
                         
    if IsZero(p) then
        return false;
    fi;

    local n, x, y, cf, A, sf;   
    x := IndeterminateOfLaurentPolynomial(p);

    # Take the largest square free divisor.
    sf := p / Gcd(p, Derivative(p));
 
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
    
    # Check if the polynomial is now 1.
    if Degree(sf) = 0 then
        return true;
    fi;

    # Check if the polynomial has even degree and is self reciprocal.
    if Degree(sf) mod 2 <> 0 or not(IsSelfReciprocalUnivariatePolynomial(sf)) then
        return false;
    fi;
    
    # Compute the number of roots in the unity sphere.
    # If this number equals the polynomial degree, then it is Kronecker.
    # Otherwise it isn't.
    y:=X(Rationals,"y");
    A:=Resultant(sf,x^2-y*x+1,x);
    # we convert A to an univariate polynomial
    A:=Value(A,[x],[1]);
    if Degree(sf) = 2*NumberOfRootsOfUnivariatePolynomialInInterval(A,-2,2) then
        return true;
    else
        return false;
    fi;
end;
