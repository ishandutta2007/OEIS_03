\\ a(n) = Sum_{k=0..floor(n/4)} binomial(2*k,2*n-8*k).
a(n) = sum(k=0, n\4, binomial(2*k, 2*n-8*k));
for(n=0, 51, print1(a(n),", "))



