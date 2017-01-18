bound:=function(p)
  local d,p21,p1;

  d:= DegreeOfUnivariateLaurentPolynomial(p);
  p21:=Value(Derivative(Derivative(p)),1);
  p1:=Value(p,1);
  return CeilingOfRational((12*p21/p1-3*d^2+6*d)/d);
end;

isKroneckerpieter:=function(p)
  local d, b,f,gd,x;
  f:=p;
  x:=IndeterminateOfLaurentPolynomial(f);
  d:=2;
  while true do
    b:=bound(f);
    if d>b then
      return false;
    fi;
    gd:=Gcd(f,x^d-1);
    if IsZero(f-gd) then
      return true;
    fi;
    if IsOne(gd) then
      d:=d+1;
    fi;
    f:=f/gd;
  od;
end;


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


isKroneckersturm:=function(p)
  local n, x, y, cf, A, psf;

  if not(IsUnivariatePolynomial(p)) then
    Error("The argument must be a polynomial in one variable");
  fi;

  if IsZero(p) then
    return false;
  fi;

  # Take the largest square free divisor.
  psf:=p/Gcd(p,Derivative(p));
  x:=IndeterminateOfLaurentPolynomial(p);

  # Remove the factors x, x+1 and x-1.
  if Value(psf, 0) = 0 then
      psf := psf / x;
  fi;
  if Value(psf, 1) = 0 then
      psf := psf/(x-1);
  fi;
  if Value(psf, -1) = 0 then
      psf := psf/(x+1);
  fi;

#  while cf[1]=0 do
#    cf:=cf{[2..Length(cf)]};
#  od;

  # Check if the polynomial is now 1.
  if Degree(psf) = 0 then
    return true;
  fi;

  # Check if the polynomial has even degree.
  if Degree(psf) mod 2 <> 0 then
    return false;
  fi;

  # Check if the polynomial is self-reciprocal.
  cf:=CoefficientsOfUnivariatePolynomial(psf);
  if not(cf=Reversed(cf)) then
    Info(InfoNumSgps,2, "The squarefreee part of the polynomial is not self-reciprocal");
    return false;
  fi;

  # Compute the number of roots in the unity sphere.
  # If this number equals the polynomial degree, then it is Kronecker.
  # Otherwise it isn't.
  y:=X(Rationals,"y");
  A:=Resultant(psf,x^2-y*x+1,x);
  # we convert A to an univariate polynomial
  A:=Value(A,[x],[1]);
  if Degree(psf) = 2*NumberOfRootsOfUnivariatePolynomialInInterval(A,-2,2) then
    return true;
  else
    return false;
  fi;
end;

isKroneckerGraeffe := function(p)
  local x, psf, p1, ps, pn, pp, coeffs_ps, factors_graeffe;

  if not(IsUnivariatePolynomial(p)) then
    Error("The argument must be a polynomial in one variable");
  fi;

  if IsZero(p) then
    return false;
  fi;

  # Take the largest square free divisor.
  psf:=p/Gcd(p,Derivative(p));
  x:=IndeterminateOfLaurentPolynomial(p);

  # Remove the factors x, x+1 and x-1.
  if Value(psf, 0) = 0 then
      psf := psf / x;
  fi;
  if Value(psf, 1) = 0 then
      psf := psf/(x-1);
  fi;
  if Value(psf, -1) = 0 then
      psf := psf/(x+1);
  fi;

  # Check if the polynomial is now 1.
  if Degree(psf) = 0 then
    return true;
  fi;

  # Check if the polynomial has even degree.
  if Degree(psf) mod 2 <> 0 then
    return false;
  fi;

  p1:= GraeffePolynomial(psf);
  if psf = p1 then
    return true;
  fi;

  if Value(psf,-x) = p1 then
    return true;
  fi;

  ps := Value(Gcd(Derivative(p1), p1), x^2);
  pp := Gcd(psf / ps, p1);
  pn := psf / (ps * pp);

#  if pp = psf then
#    return true;
#  fi;

  factors_graeffe := Difference([ps, pp, pn], [1+0*x]);

  if Length(factors_graeffe) = 1 then
    if IsOne(ps) then
      return false;
    else
      return isKroneckerGraeffe(p1 /  Gcd(p1,Derivative(p1)));
    fi;
  else
    return ForAll(factors_graeffe, isKroneckerGraeffe);
  fi;
end;
