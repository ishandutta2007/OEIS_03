seq(n) = my(A=1); for(i=1, n, A=1 + x/A^3 * (1 + 1/A^3) +x*O(x^n)); Vec(A);      
seq(35)    

seq(n) = my(A=1); for(i=1, n, A=1 + x/A^3 * (1 + 1/A^3) +x*O(x^n)); Vec(subst(1/A, x, -x));
seq(35)      