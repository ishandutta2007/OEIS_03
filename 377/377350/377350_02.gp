\\ E.g.f. satisfies A(x) = 1 - log(1 - x*A(x)^3)/A(x)^3.
my(A=1, n=32); for(i=1, n, A=1 - log(1 - x*A^3)/A^3 +x*O(x^n) ); Vec(serlaplace(A))


