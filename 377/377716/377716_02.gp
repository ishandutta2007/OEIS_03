\\ E.g.f. satisfies A(x) = (1 + (exp(x) - 1) * A(x))^2.
my(A=1, n=25); for(i=1, n, A=(1 + (exp(x + x*O(x^n)) - 1) * A)^2 ); Vec(serlaplace(A))

my(A=1, n=25); for(i=1, n, A=(1 + (exp(x + x*O(x^n)) - 1) * A)^2 ); Vec(serlaplace(A^(1/2)))

