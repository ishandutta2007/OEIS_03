M=18;

a(n) = sum(k=1, n, sumdiv(k, d, kronecker(-4, k/d)*(d^n-(d-1)^n)));
for(n=1, M, print1(a(n), ", "));
