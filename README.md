# MapleSturmDiscriminants
A Maple package for the computation of Sturm discriminants


The follwoing has been implemented:

1. Create a function that generates a degree d Sturm generic
polynomial.

2. Check if the coefficients of the elimination ideals are
algebraically independent, for example, using the Jacobian Criterion.

3. If 2 is true, then susbstitute the computation of the Sturm
sequence of the elimination ideal by the computation of the Sturm
sequence of a generic polynomial of the same degree.

4. If 2 is true, substitute back the coefficients of the elimination
ideal.


Apparently, step 4 is the slowest one. I suspect that this is due to
automatic simplifications in Maple. I will write a function that
avoids such automatic simplifications