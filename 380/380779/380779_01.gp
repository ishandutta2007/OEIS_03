M=17;

a(n, q=1, r=1, s=1, t=-1, u=2) = q*n!*sum(k=0, n, (r*n+(s-r)*k+q)^(k-1)*binomial(r*u*n+((s-r)*u+t)*k+q*u, n-k)/k!);
for(n=0, M, print1(a(n),", "));

\\ a(n) = n! * Sum_{k=0..n} (n+1)^(k-1) * binomial(2*n-k+2,n-k)/k!.
b(n) = n!*sum(k=0, n, (n+1)^(k-1)*binomial(2*n-k+2,n-k)/k!);
for(n=0, 100, print1(a(n)-b(n),", "));


