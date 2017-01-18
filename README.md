# KroneckerPolynomials

Algorithms to check whether a polynomial is Kronecker or not. The implementation is carried out in GAP, and will susbstitute the current implementation of IsKroneckerPolynomial in numericalsgps [2](#numericalsgps). You can find the following algorithms:

1. Sequential GCD method. Credits to Pieter Moree.
2. Sturm sequence based method. Credits to David Boyd.
3. Graeffe polynomial based method. See [1](#davenport), page 248.

## References

<a name="davenport">1</a>: Effective Tests for Cyclotomic Polynomials, R. J. Bradford and J.H. Davenport, 1988, International Symposium on Symbolic and Algebraic Computation (pp. 244-251). Springer Berlin Heidelberg.

<a name="numericalsgps">2</a>: Delgado, M., Garcia-Sanchez, P. A. and Morais, J., NumericalSgps, A
package  for  numerical  semigroups, Version 1.0.1 dev (2015), (Refereed GAP
package), http://www.fc.up.pt/cmup/mdelgado/numericalsgps/.
