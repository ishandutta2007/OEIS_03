\\ G.f. A(x) satisfies A(x) = ( 1 + x/(1 - x*A(x)^(2/3)) )^3.
seq(n) = my(A=1); for(i=1, n, A=(1 + x/(1 - x*A^(2/3) +x*O(x^n)))^3 ); Vec(A);      
seq(18)

\\ G.f.: A(x) = (1 + x*B(x))^3 where B(x) is the g.f. of A161634.
seq(n) = my(A=1); for(i=1, n, A=(1 + x/(1 - x*A^(2/3) +x*O(x^n)))^3 ); Vec(A^(1/3));
seq(18)
