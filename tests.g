pieter:=function(s)
    local g,gaps;

    gaps:=GapsOfNumericalSemigroup(s);
    g:=GenusOfNumericalSemigroup(s);
    return 1+2*g-4*Length(Filtered(gaps, x->x mod 2=0));



end;

pieter2:=function(s)
    local godd,geven,gaps;

    gaps:=GapsOfNumericalSemigroup(s);
    godd:=Length(Filtered(gaps,x->x mod 2<>0));
    geven:=Length(Filtered(gaps,x-> x mod 2=0));
    return godd<geven;


end;

pietersearch:=function(n)
    local lc,li,ld;
    li:=IrreducibleNumericalSemigroupsWithFrobeniusNumber(n);
    lc:=CompleteIntersectionNumericalSemigroupsWithFrobeniusNumber(n);
    ld:=Difference(li,lc);
    return Filtered(ld, s->not(pieter2(s)));

end;

pietersearchupto:=function(n)
    local lc,lct,li,lit,ld,lp,diff, i;
    lc:=[];
    li:=[];
    ld:=[];
    lp:=[];
    for i in [1..n] do
        lit:=IrreducibleNumericalSemigroupsWithFrobeniusNumber(2*i-1);
        li:=Union(li,lit);
        lct:=CompleteIntersectionNumericalSemigroupsWithFrobeniusNumber(2*i-1);
        lc:=Union(lc,lct);
        diff:=Difference(lit,lct);
        ld:=Union(ld,diff);
        lp:=Union(lp,Filtered(diff, s->not(pieter2(s))));
        if (Length(ld)>0) then
            Print(2*i-1," -> ",Float(Length(lp)/Length(ld)),"\n");
        fi;
    od;
    return Length(lp)/Length(ld);
end;



pieteraverage:=function(n)
    local lc,li,ld;
    li:=IrreducibleNumericalSemigroupsWithFrobeniusNumber(n);
    lc:=CompleteIntersectionNumericalSemigroupsWithFrobeniusNumber(n);
    ld:=Difference(li,lc);
    return Length(Filtered(ld, s->not(pieter2(s))))/Length(ld);

end;

sem_test:=function(m,q)
	local gen;

	if m<2*q+3 then
		Error ("Check m and q");
	fi;

	gen:=Union([m,m+1],[(q*m+2*q+2)..(q*m+m-1)]);
	return NumericalSemigroup(gen);
end;


family:=function(m)
	if (m-3) mod 2 =0 then
		return List([1..(m-3)/2], x->sem_test(m,x));
	else
		return List([1..(m-4)/2],  x->sem_test(m,x));
	fi;
end;


diameter:=function(s)
	local gs;
	gs:=GapsOfNumericalSemigroup(s);
	return(Maximum(Set([2..GenusOfNumericalSemigroup(s)], i ->gs[i]-gs[i-1])));
end;

denominator:=function(s,x)
	local msg;
	msg:=MinimalGeneratingSystemOfNumericalSemigroup(s);
	return SemigroupPolynomial(s,x)*Product(List(msg, n->(1-x^n)));
end;

conjecture:=function(s,x)
	local bet;
	bet:=BettiElementsOfNumericalSemigroup(s);
	return denominator(s,x)=(1-x)*Product(List(bet, n->(1-x^n)^(Length(FactorizationsElementWRTNumericalSemigroup(n,s))-1)));
end;

newden:=function(x, d,s)
	local i, out, el, k;

	i:=0;
	el:=x;
	k:=DenumerantElementInNumericalSemigroup(d,s)-1;
	out :=DenumerantElementInNumericalSemigroup(el,s);

	while (el-d) in s do
		el:=el-d;
		i:=i+1;
		out:=out+(-1)^i*Binomial(k,i)*DenumerantElementInNumericalSemigroup(el,s);
	od;
	return out;
end;

sumshaded:=function(x,s)
	return EulerCharacteristic(ShadedSetOfElementInNumericalSemigroup(x,s));
end;

denominator2:=function(s,x)
	local msg, f, bound, l;

	msg:=MinimalGeneratingSystemOfNumericalSemigroup(s);
	f:=FrobeniusNumber(s);

	bound:=f+Sum(msg);
	l:=FirstElementsOfNumericalSemigroup(bound,s);
	return (1-x)*Sum(List(l, n-> sumshaded(n,s)*x^n));
end;

fk:=function(k,x)
  return 1-x+x^k-x^(2*k-1)+x^(2*k);
end;
