#############################################################################
#
#  test_kronecke_fast_3.sh        Andrés Herrera-Poyatos <andreshp9@gmail.com>
#                                Pedro A. Garcia-Sanchez <pedro@ugr.es>
#
# Test the kronecker algorithms with the polynomial \Phi_n(x).
# Only the fastest algorihtms are tested.
#
#############################################################################

gap -r -b -q << EOI
LoadPackage("num");;
Read("kronecker_pieter.g");;
Read("kronecker_sturm.g");;
Read("kronecker_graeffe.g");;

x := X(Rationals, "x");;
m := 250;;

Print("n,Pieter,Sturm,Graeffe\n");

for i in [1..m] do
    Print(i);
    Print(",");
    p := CyclotomicPolynomial(Rationals, i);
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
    if not(sol_pieter) then
        Print("Error in current algorithm, the polynomial ", p, " is Kronecker.");
    fi;
    if not(sol_sturm) then
        Print("Error in current algorithm, the polynomial ", p, " is Kronecker.");
    fi;
    if not(sol_graeffe) then
        Print("Error in current algorithm, the polynomial ", p, " is Kronecker.");
    fi;
od;

quit;
EOI
