\\ G.f. A(x) satisfies A(x) = 1/( 1 - x*A(x)/(1 - x*A(x)) )^6.
my(A=1, n=22); for(i=1, n, A= 1/(1 - x*A/(1 - x*A))^6 + x*O(x^n) ); Vec(A)

\\ G.f.: A(x) = B(x)^6 where B(x) is the g.f. of A378694.
my(A=1, n=22); for(i=1, n, A= 1/(1 - x*A/(1 - x*A))^6 + x*O(x^n) ); Vec(A^(1/6))

 