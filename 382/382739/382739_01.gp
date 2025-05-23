M=15;

\\ A(n,k) = Sum_{j=0..min(n,k)} (j!)^2 * binomial(j+3,3) * Stirling2(n,j) * Stirling2(k,j).
a(n, k) = sum(j=0, min(n,k), (j!)^2 * binomial(j+3,3) * stirling(n,j,2) * stirling(k,j,2));

\\ a(n) = Sum_{k=0..n} (k!)^2 * binomial(k+3,3) * Stirling2(n,k)^2.
b(n) = sum(k=0, n, (k!)^2 * binomial(k+3,3) * stirling(n,k,2)^2);
for(n=0, M, print1(b(n), ", "));

\\ a(n) == 0 (mod 4) for n > 0.
for(n=0, 100, print1(b(n)%4, ", "));

\\ Main diagonal of A382736.
for(n=0, 100, print1(a(n, n)-b(n), ", "));




