M=15;

a(n) = polcoef(1/prod(k=1, 2*n, 1-k*x+x*O(x^n)), n);
for(n=0, M, print1(a(n), ", "));