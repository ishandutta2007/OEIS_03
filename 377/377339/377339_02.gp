\\ E.g.f. satisfies A(x) = ( 1 + (exp(x*A(x)) - 1)/A(x) )^2.
my(A=1, n=32); for(i=1, n, A=(1 + (exp(x*A) - 1)/A)^2 +x*O(x^n) ); Vec(serlaplace(A))

my(A=1, n=32); for(i=1, n, A=(1 + (exp(x*A) - 1)/A)^2 +x*O(x^n) ); Vec(serlaplace(A^(1/2)))

