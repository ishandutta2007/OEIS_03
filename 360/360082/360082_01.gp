M=22;

a(n) = sum(k=0, n, binomial(4*k, n-k)*binomial(2*k, k)/(k+1));
for(n=0, M, print1(a(n), ", "));