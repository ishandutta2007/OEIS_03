seq(n) = my(A=1); for(i=1, n, A=(1 + x * A^(1/5) * (1+A) +x*O(x^n) )^(5/2) ); Vec(A);
seq(25)    

\\ A^(1/5) = A370472
seq(n) = my(A=1); for(i=1, n, A=(1 + x * A^(1/5) * (1+A) +x*O(x^n) )^(5/2) ); Vec(A^(1/5));
seq(25)
