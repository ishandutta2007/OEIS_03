M=21;

a(n) = sum(k=0, n\3, (3*k)!*abs(stirling(n, 3*k, 1))/k!);
for(n=0, M, print1(a(n), ", "));
