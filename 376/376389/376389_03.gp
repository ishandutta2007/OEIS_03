M=17;

\\ E.g.f. A(x) satisfies A(x) = 1/(2 - exp(x*A(x)))^2.
my(A=1, n=M); for(i=1, n, A=1/(2 - exp(x*A + x*O(x^n)))^2); Vec(serlaplace(A))

\\ E.g.f.: B(x)^2, where B(x) is the e.g.f. of A367134.
my(A=1, n=M); for(i=1, n, A=1/(2 - exp(x*A + x*O(x^n)))^2); Vec(serlaplace(A^(1/2)))
