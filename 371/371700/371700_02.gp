seq(n) = my(A=1); for(i=1, n, A=1 + x * A^6 * (1+A) +x*O(x^n) ); Vec(A);
seq(25)    
