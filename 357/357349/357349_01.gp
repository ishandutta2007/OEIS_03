M=18;

a(n) = sum(k=1, n, (n+k)^(k-1)*stirling(n, k, 1));
for(n=0, M, print1(a(n), ", "));