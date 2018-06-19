#############################################################################
#
#  test_kronecke_fast_4.sh        Andrés Herrera-Poyatos <andreshp9@gmail.com>
#                                Pedro A. Garcia-Sanchez <pedro@ugr.es>
#
# Test the kronecker algorithms with the polynomial (x^n+2)(x^n-1).
# Only the fastest algorihtms are tested.
#
#############################################################################

gap.sh -r -b -q << EOI
LoadPackage("num");;
Read("kronecker_pieter.g");;
Read("kronecker_sturm.g");;
Read("kronecker_graeffe.g");;

x := X(Rationals, "x");;
m := 16;;

Print("n,Pieter (improved),Sturm,Graeffe (improved)\n");

for i in [1..m] do
    Print(i);
    Print(",");
    p := (x^(2^i)+2)*(x^(2^i)-1);
    start := Runtime();
    sol_pieter := IsKroneckerPolynomialPieterImproved(p);
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
    if sol_pieter then
        Print("Error in Pieter algorithm, the polynomial ", p, " is not Kronecker.");
    fi;
    if sol_sturm then
        Print("Error in Sturm algorithm, the polynomial ", p, " is not Kronecker.");
    fi;
    if sol_graeffe then
        Print("Error in Graeffe algorithm, the polynomial ", p, " is not Kronecker.");
    fi;
od;

quit;
EOI
