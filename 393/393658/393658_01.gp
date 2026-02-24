M=33;

\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/6)} binomial(n+k,k) * binomial(n+1,n-6*k).
a(n) = (1/(n+1)) * sum(k=0, n\6, binomial(n+k,k)*binomial(n+1,n-6*k));
for(n=0, M, print1(a(n),", "))  



