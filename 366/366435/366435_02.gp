seq(n) = my(A=1); for(i=1, n, A=1+x * ((1 - x +x*O(x^n))/A)^(5/2) ); Vec(A); 
seq(30)