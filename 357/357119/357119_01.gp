M=11;

T(n, k) = sum(j=0, n, abs(stirling(n, k*j, 1)));
for(n=0, M, for(k=0, n, print1(T(k, n-k), ", ")));