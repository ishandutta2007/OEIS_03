\\ a(n) = (n+1) * (n!)^3 * Sum_{k=1..n} 1/((k+1) * (k!)^3).
a(n) = (n+1)*n!^3*sum(k=1, n, 1/((k+1)*k!^3));
for(n=0, 13, print1(a(n), ", "));

\\ a(n) = A368775(n) - (n+1) * (n!)^3.
A368775(n) = (n+1)*n!^3*sum(k=0, n, 1/((k+1)*k!^3));
for(n=0, 13, print1(A368775(n), ", "));
a(n) = A368775(n) - (n+1)*n!^3;
for(n=0, 16, print1(a(n), ", "));
