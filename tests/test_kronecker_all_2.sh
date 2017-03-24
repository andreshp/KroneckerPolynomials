#############################################################################
##
#W  test_kronecke_all_2.sh       Andrés Herrera-Poyatos <andreshp9@gmail.com>
#W                               Pedro A. Garcia-Sanchez <pedro@ugr.es>
##
## Test the kronecker algorithms.
##
#############################################################################

gap -r -b -q << EOI
LoadPackage("num");;
Read("kronecker_pieter.g");;
Read("kronecker_sturm.g");;
Read("kronecker_graeffe.g");;

x := X(Rationals, "x");;
m := 9;;

Print("n,NumerisalSgp,Pieter,Sturm,Graeffe\n");
    
for i in [1..m] do
    Print(i);
    Print(",");
    p := x^(2^i)-2;
    start := Runtime();
    sol_current := IsKroneckerPolynomial(p);
    total:=Runtime()-start;
    Print(total);
    Print(",");
    start := Runtime();
    sol_pieter := IsKroneckerPolynomialPieter(p);
    total:=Runtime()-start;
    Print(total);
    Print(",");
    start := Runtime();
    sol_sturm := IsKroneckerPolynomialSturm(p);
    total:=Runtime()-start;
    Print(total);
    Print(",");
    start := Runtime();
    sol_graeffe := IsKroneckerPolynomialGraeffe(p);
    total:=Runtime()-start;
    Print(total);
    Print("\n");
    if sol_current then
        Print("Error in current algorithm, the polynomial ", p, " is not Kronecker.");
    fi;
    if sol_pieter then
        Print("Error in current algorithm, the polynomial ", p, " is not Kronecker.");
    fi;
    if sol_sturm then
        Print("Error in current algorithm, the polynomial ", p, " is not Kronecker.");
    fi;
    if sol_graeffe then
        Print("Error in current algorithm, the polynomial ", p, " is not Kronecker.");
    fi;
od;

quit;
EOI
