M=20;

a(n) = sum(k=0, n, (-1)^k*(3*k)^(n-k)*binomial(n, k));
for(n=0, M, print1(a(n), ", "));