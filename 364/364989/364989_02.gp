seq(n) = my(A=1); for(i=1, n, A=1+x*A^4*exp(x*A^4 +x*O(x^n)) ); Vec(serlaplace(A));         
seq(18)           

