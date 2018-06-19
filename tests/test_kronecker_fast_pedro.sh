#############################################################################
#
#  test_kronecke_all_pedro.sh        Andrés Herrera-Poyatos <andreshp9@gmail.com>
#                                Pedro A. Garcia-Sanchez <pedro@ugr.es>
#
# Test the kronecker algorithms.
#
#############################################################################

gap.sh -r -b -q << EOI
LoadPackage("num");;
Read("kronecker_pieter.g");;
Read("kronecker_sturm.g");;
Read("kronecker_graeffe.g");;

x := X(Rationals, "x");;
m := 200;;

Print("n,Pieter (improved),Boyd,Graeffe (improved)\n");

for i in [1..m] do
    Print(i);
    Print(",");
    p := 1 - x + x^i - x^(2*i-1) + x^(2*i);
    start := Runtime();
    sol_pieter := IsKroneckerPolynomialPieterImproved(p);
    total:=Runtime()-start;
    Print(total);
    Print(",");
    start := Runtime();
    sol_pieter_imp := IsKroneckerPolynomialSturm(p);
    total:=Runtime()-start;
    Print(total);
    Print(",");
    start := Runtime();
    sol_graeffe := IsKroneckerPolynomialGraeffe(p);
    total:=Runtime()-start;
    Print(total);
    Print("\n");
od;

quit;
EOI
