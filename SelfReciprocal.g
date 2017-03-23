IsSelfReciprocal := function(p)
    local cf;
    
    if not(IsUnivariatePolynomial(p)) then
	Error("The argument must be a polynomial\n");
    fi;
    
    # Check if the polynomial is self-reciprocal.
    cf:=CoefficientsOfUnivariatePolynomial(p);
    return cf=Reversed(cf);    
end;


IsSelfReciprocalTest := function()
    local y;
    y := X(Rationals, "y");
    for i in [2..1000] do
        IsSelfReciprocal(y^i+1);
        IsSelfReciprocal(y^i-1);
        IsSelfReciprocalUnivariatePolynomial(y^i+y^Int(i/2)+1);
        IsSelfReciprocalUnivariatePolynomial(y^i+y^Int(i/2)+y^(Int(i/2)+1)+1);
    od;
end;

IsSelfReciprocalUTest := function()
    local y;
    y := X(Rationals, "y");
    for i in [2..1000] do
        IsSelfReciprocalUnivariatePolynomial(y^i+1);
        IsSelfReciprocalUnivariatePolynomial(y^i-1);
        IsSelfReciprocalUnivariatePolynomial(y^i+y^Int(i/2)+1);
        IsSelfReciprocalUnivariatePolynomial(y^i+y^Int(i/2)+y^(Int(i/2)+1)+1);
    od;
end;
