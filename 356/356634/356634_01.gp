M=22;

a(n) = n!*sum(k=0, n\4, (n-4*k)^k/24^k);
for(n=0, M, print1(a(n), ", "));