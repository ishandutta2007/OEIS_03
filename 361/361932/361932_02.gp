seq(n) = my(A=1); for(i=1, n, A=1 + x^3*A^3/(1-x +x*O(x^n))^3 ); Vec(A); 
seq(30)
         