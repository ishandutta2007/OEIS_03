M=30;

a(n) = sum(k=0, n, k!*binomial(n, k)*binomial(2*n+1, k));
for(n=0, M, print1(a(n), ", "));
