M=15;

a(n, q=3, r=3, s=3, t=0, u=1) = q*n!*sum(k=0, n, (r*n+(s-r)*k+q)^(k-1)*binomial((r*u+1)*n+((s-r)*u+t-1)*k+q*u-1, n-k)/k!);
for(n=0, M, print1(a(n),", "));

\\ a(n) = 3 * n! * Sum_{k=0..n} (3*n+3)^(k-1) * binomial(4*n-k+2,n-k)/k!.
b(n) = 3 * n! * sum(k=0, n, (3*n+3)^(k-1) * binomial(4*n-k+2,n-k)/k!);
for(n=0, 50, print1(a(n)-b(n),", "));

