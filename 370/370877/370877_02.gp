M=19;

\\ a(n) = n! * Sum_{k=0..floor(n/2)} (2*k+1)^(k-1) * binomial(n,2*k)/(2^k * k!).
a(n) = n! * sum(k=0, n\2, (2*k+1)^(k-1) * binomial(n,2*k)/(2^k * k!));
for(n=0, M, print1(a(n), ", "))