M=20;

a(n) = n!*sum(k=0, n\5, abs(stirling(n-4*k, k, 1))/(n-4*k)!);
for(n=0, M, print1(a(n), ", "));
