\\ E.g.f. satisfies A(x) = (1 + x * A(x) * exp(x*A(x)))^2.
my(A=1, n=25); for(i=1, n, A=(1 + x * A * exp(x*A +x*O(x^n)))^2 ); Vec(serlaplace(A))

my(A=1, n=25); for(i=1, n, A=(1 + x * A * exp(x*A +x*O(x^n)))^2 ); Vec(serlaplace(A^(1/2)))
