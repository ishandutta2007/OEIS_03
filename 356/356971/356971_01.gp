M=20;

a(n) = n!*sum(k=0, n\4, (n-2*k+1)^(k-1)*abs(stirling(n-3*k, k, 1))/(n-3*k)!);
for(n=0, M, print1(a(n), ", "));