\\ a(n) = (1/(n+1)) * Sum_{k=0..floor(n/2)} binomial(n+k,k) * binomial(5*n-2*k+3,n-2*k).
a(n) = (1/(n+1)) * sum(k=0, n, binomial(n+k,k) * binomial(5*n-2*k+3,n-2*k) );
for(n=0, 22, print1(a(n),", "))  

