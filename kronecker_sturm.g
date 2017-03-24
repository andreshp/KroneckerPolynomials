#############################################################################
##
#W  kronecker_sturm.g           Andrés Herrera-Poyatos <andreshp9@gmail.com>
#W                              Pedro A. Garcia-Sanchez <pedro@ugr.es>
##
## Algorithm to check whether a polynomial is Kronecker based on
## the Sturm sequence and Sturm's Theorem.
##
#############################################################################

########################################################################
##
#F SturmSequence(f) returns the Sturm sequence of the polynomial f.
##
#########################################################################
#DeclareGlobalFunction("SturmSequence");

########################################################################
##
#F NumberOfRootsOfUnivariatePolynomialInInterval(f, a, b) returns the 
## number of roots that the polynomial f has in the interval [a,b].
## It uses the sturm sequence of the polynomial.
##
#########################################################################
#DeclareGlobalFunction("NumberOfRootsOfUnivariatePolynomialInInterval");

########################################################################
##
#F IsKroneckerPolynomialSturm(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials
##
#########################################################################
#DeclareGlobalFunction("IsKroneckerPolynomialSturm");

########################################################################
##
#F SturmSequence(f) returns the Sturm sequence of the polynomial f.
##
#########################################################################
SturmSequence := function(f)
    local sst;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    sst := [f, Derivative(f)];
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
NumberOfRootsOfUnivariatePolynomialInInterval := function(f,a,b)
    local sst, sgnc, valuesa, valuesb, sf;

    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    # Take the largest square free divisor.
    sf:=f/Gcd(f,Derivative(f));
    
    # Counts the number of sign changes in the sequence of numbers given
    # as argument.
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
    
    # Apply Sturm theorem to obtain the number of zeroes in [a,b].
    sst := SturmSequence(sf);
    valuesa := List(sst, p -> Value(p, a));
    valuesb := List(sst, p -> Value(p, b));

    return sgnc(valuesa) - sgnc(valuesb);
end;

########################################################################
##
#F IsKroneckerPolynomialSturm(f) decides if
##   f is a Kronecker polynomial, that is, a monic polynomial with integer 
##   coefficients having all its roots in the unit circunference, equivalently, 
##   is a product of cyclotomic polynomials.
##
##   The algorithm was proposed by David Boyd. It uses the sturm's algorithm 
##   to count roots.
##
##   First we obtain the polynomial sf as the largest square free factor of f.
##   Then we check whether x, x-1 or x+1 divide sf. If it is the case,
##   then we remove these factors from sf. Afterwards, let n be the degree of sf. 
##   If n is 0, then  f is Kronecker. If n is not even, then f is not Kronecker.
##   If sf is not Self-Reciprocal then f is not Kronecker neither.
##   Now we can write uniquely sf(x) = x^(n/2) A(x+1/x), where A is a polynomial
##   of degree n/2. A can be computed as the resultant of sf and x^2-y*x+1.
##   Note that A has half as many roots in [-2,2] as sf has roots in the unit 
##   circle. Now we use the sturm sequence algorithm to determine the 
##   number of roots of A in [-2,2]. f is Kronecker if and only if this
##   number equals n.
#########################################################################
IsKroneckerPolynomialSturm := function(f)
    local x, y, A, sf;
    if not(IsUnivariatePolynomial(f)) then
        Error("The argument must be a polynomial in one variable");
    fi;

    if IsZero(f) then
        return false;
    fi;

    x := IndeterminateOfLaurentPolynomial(f);

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

    # Compute the number of roots in the unit circunference.
    # If this number equals the polynomial degree, then it is Kronecker.
    # Otherwise it isn't.
    y:=X(Rationals, "y");
    A:=Resultant(sf, x^2-y*x+1, x);
    # we convert A to an univariate polynomial
    A:=Value(A, [x], [1]);

    return Degree(sf) = 2*NumberOfRootsOfUnivariatePolynomialInInterval(A,-2,2);
end;
