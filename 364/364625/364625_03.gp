seq(n) = my(A=1); for(i=1, n, A= 1/(1 - x +x*O(x^n))^3 + x^2*A^2 ); Vec(A); 
seq(20)                                                                             
