M=24;

a(n) = sum(k=0, n\3, 2^k*stirling(n, 3*k, 2));
for(n=0, M, print1(a(n), ", "));