seq(n) = my(A=1); for(i=1, n, A=1 + x*A^4/(1 - x*A +x*O(x^n)) ); Vec(A);       
seq(30)   

\\ G.f. A(x) satisfies A(x) = 1/(1 - x*A(x)^3/(1 - x*A(x))).
seq(n) = my(A=1); for(i=1, n, A=1/(1 - x*A^3/(1 - x*A) +x*O(x^n)) ); Vec(A);
seq(30)
