seq(n) = my(A=1); for(i=1, n, A=1 + x*A*(1 + x^4*A^5 +x*O(x^n)) ); Vec(A); 
seq(40)  

