M=24;

\\ a(n) = n! * Sum_{k=0..floor(n/4)} binomial(n/4-k,k)/(n-4*k)!.
a(n) = n! * sum(k=0, n\4, binomial(n/4-k,k)/(n-4*k)!);
for(n=0, M, print1(a(n), ", "))