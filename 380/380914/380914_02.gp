\\ E.g.f. A(x) satisfies A(x) = exp(x / (1 - x*A(x))) / (1 - x*A(x)).
my(A=1, n=22); for(i=1, n, A=exp(x / (1 - x*A)) / (1 - x*A) + x*O(x^n) ); Vec(serlaplace(A))
