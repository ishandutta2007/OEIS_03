\\ E.g.f. A(x) satisfies A(x) = exp(2*x*A(x)) / ( 1 - x*exp(x*A(x)) ).
my(A=1, n=25); for(i=1, n, A = exp(2*x*A) /(1 - x * exp(x * A + x*O(x^n))) ); Vec(serlaplace(A))

