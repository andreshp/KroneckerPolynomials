#########################################
# SturmSequence(p)
# returns the sturm sequence of the polynomial p
#########################################
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

#########################################
# SemigroupPolynomial(s,x)
# Counts the number of roots of an univariate polynomial in the given interval.
# It uses the sturm sequence.
#########################################
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
