M=25;

\\ a(n) = Sum_{k=0..floor(n/3)} binomial(n+k,k)^2 * binomial(n-2*k,k).
a(n) = sum(k=0, n\3, binomial(n+k,k)^2 * binomial(n-2*k,k));
for(n=0, M, print1(a(n),", "))         