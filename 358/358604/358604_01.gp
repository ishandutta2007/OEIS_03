M=37;

a(n) = sum(k=0, n\3, (-1)^k*(n-2*k)!/(n-3*k)!);
for(n=0, M, print1(a(n), ", "));