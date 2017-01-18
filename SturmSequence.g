SturmSequence := function(p)
    local sst;
    sst := [p, Derivative(p)];
    while Degree(sst[Length(sst)]) >= 1 do
        Add(sst, -QuotientRemainder(sst[Length(sst)-1],sst[Length(sst)])[2]);
    od;
    return sst;
end;

NumberOfRootsOfUnivariatePolynomialInInterval := function(p,a,b)
    local sst, sgnc, valuesa, valuesb;

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
    valuesa := List(sst, p -> Value(p, a));
    valuesb := List(sst, p -> Value(p, b));
    Print(valuesa);
    Print(valuesb);

    return sgnc(valuesa) - sgnc(valuesb);
end;
