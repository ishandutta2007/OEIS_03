a(n, k, m) = my(A=1, B=sum(k=0, n, binomial(3*k,k)/(2*k+1) * x^k) + x*O(x^n)); for(i=1, n, A=1/(1 - x*A^k) * subst(B, x, x*A^m) ); Vec(A);
a(20, 1, 1)
