M=11;

T(n, k) = sum(j=1, n, sumdiv(j, d, (-1)^(j/d+1)*(d^k-(d-1)^k)));
for(n=1, M, for(k=1, n, print1(T(k, n-k+1), ", ")));
