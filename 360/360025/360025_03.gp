M=20;

a(n) = if(n==0, 1, 1 - sum(k=0, n-3, a(k) * a(n-k-3)));
for(n=0, M, print1(a(n), ", "));