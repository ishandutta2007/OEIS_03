seq(n) = my(A=1); for(i=1, n, A=1/(1 - x*A +x*O(x^n))^A^3 ); Vec(serlaplace(A));         
seq(30)    


