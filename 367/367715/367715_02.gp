seq(n) = my(A=1); for(i=1, n, A=1/((1 - x) * (1 - x * (1 + x + x^2 + x^3) * subst(A, x, x^4) +x*O(x^n)))); Vec(A); 
seq(30)        