\\ a(n) = (1/n) * Sum_{k=0..floor((n-1)/2)} 4^(n-k) * binomial(n,k) * binomial(3*n-k,n-1-2*k) for n > 0.
a(n) = if(n==0, 1, sum(k=0, (n-1)\2, 4^(n-k)*binomial(n, k)*binomial(3*n-k, n-1-2*k))/n);
for(n=0, 17, print1(a(n),", ")) 

\\ a(n) = 2^n * A219538(n).
for(n=0, 20, print1(a(n)/2^n, ", "))


