seq(n) = my(A=1); for(i=1, n, A=1 + x*A^2 + x^2*A^6 +x*O(x^n)); Vec(A);  
seq(35)  
