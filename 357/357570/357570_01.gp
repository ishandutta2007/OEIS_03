M=46;

a(n) = sum(k=0, n\5, (n-4*k)!/(n-5*k)!);
for(n=0, M, print1(a(n), ", "));