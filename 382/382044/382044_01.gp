M=17;

\\ a(n) = n! * Sum_{k=0..n} (2*k)^(n-k) * binomial(n+3*k+1,k)/((n+3*k+1) * (n-k)!).
a(n) = n! * sum(k=0, n, (2*k)^(n-k) * binomial(n+3*k+1,k)/((n+3*k+1) * (n-k)!));
for(n=0, M, print1(a(n),", "));

