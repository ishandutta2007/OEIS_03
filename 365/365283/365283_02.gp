seq(n) = my(A=1); for(i=1, n, A=1 + x*A*exp(x^2*A^2 +x*O(x^n)) ); Vec(serlaplace(A));         
seq(30)        
