a(n) = (-1)^n * sum(k=0, n, 9^k * binomial(-1/3,k) * binomial(2*k,n-k));
for(n=0, 22, print1(a(n),", "))
