M=19;

\\ E.g.f. A(x) satisfies A(x) = 1/(1 - x*A(x) * (exp(x*A(x)) - 1))^2.
my(A=1, n=M); for(i=1, n, A=1/(1 - x*A*(exp(x*A + x*O(x^n)) - 1))^2); Vec(serlaplace(A))

\\ E.g.f.: B(x)^2, where B(x) is the e.g.f. of A371271.
my(A=1, n=M); for(i=1, n, A=1/(1 - x*A*(exp(x*A + x*O(x^n)) - 1))^2); Vec(serlaplace(A^(1/2)))


