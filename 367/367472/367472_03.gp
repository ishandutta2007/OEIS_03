\\ a(0) = 1; a(n) = 3*Sum_{k=1..n} (k/n + 1) * binomial(n,k) * a(n-k).
a(n) = if(n==0, 1, 3*sum(k=1, n, (k/n + 1)*binomial(n,k)*a(n-k)));
for(n=0, 17, print1(a(n),", "))   
