\\ a(n) = Sum_{k=0..n-1} binomial(5*k,k) * binomial(5*n-5*k,n-k-1).
a(n) = sum(k=0, n-1, binomial(5*k, k) * binomial(5*n-5*k, n-k-1));
for(n=0, 21, print1(a(n), ", "));


