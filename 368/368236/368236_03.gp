\\ a(0) = 1; a(n) = 2*n*a(n-1) + Sum_{k=1..n} (-1)^(k-1) * binomial(n,k) * a(n-k).
a(n) = if(n==0, 1, 2*n*a(n-1) + sum(k=1, n, (-1)^(k-1)*binomial(n, k)*a(n-k)));
for(n=0, 15, print1(a(n), ", "));

