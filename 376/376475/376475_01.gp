M=21;

\\ a(n) = n! * Sum_{k=0..floor(n/3)} (n+1)^(k-1) * binomial(n-2*k-1,n-3*k)/k!.
a(n) = n! * sum(k=0, n\3, (n+1)^(k-1) * binomial(n-2*k-1,n-3*k)/k!);
for(n=0, M, print1(a(n), ", "));