M=16;

a(n, q=2, r=2, s=2, t=0, u=1) = q*n!*sum(k=0, n, (r*n+(s-r)*k+q)^(k-1)*binomial((r*u+1)*n+((s-r)*u+t-1)*k+q*u-1, n-k)/k!);
for(n=0, M, print1(a(n),", "));

\\ a(n) = 2 * n! * Sum_{k=0..n} (2*n+2)^(k-1) * binomial(3*n-k+1,n-k)/k!.
b(n) = 2 * n! * sum(k=0, n, (2*n+2)^(k-1) * binomial(3*n-k+1,n-k)/k!);
for(n=0, 50, print1(a(n)-b(n),", "));

