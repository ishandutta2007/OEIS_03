\\ a(n) = (1/(n+1)) * Sum_{k=0..n} 5^k * (-1)^(n-k) * (1/2)^(n-2*k) * binomial(n+k,k) * binomial(k,n-k).
a(n) = sum(k=0, n, 5^k * (-1)^(n-k) * (1/2)^(n-2*k) * binomial(n+k,k) * binomial(k,n-k)) / (n+1);
for(n=0, 18, print1(a(n), ", "));
