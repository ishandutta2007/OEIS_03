M=23;

a(n, r=3, s=2, t=2, u=0) = r*sum(k=0, n, binomial(t*k+u*(n-k)+r, k)*binomial(n+(s-1)*k-1, n-k)/(t*k+u*(n-k)+r));
for(n=0, M, print1(a(n), ", "))

\\ a(n) = 3 * Sum_{k=0..n} binomial(2*k+3,k) * binomial(n+k-1,n-k)/(2*k+3).
b(n) = 3 * sum(k=0, n, binomial(2*k+3, k)*binomial(n+k-1, n-k)/(2*k+3));
for(n=0, 25, print1(a(n)-b(n), ", "))
