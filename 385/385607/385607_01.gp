M=17;

\\ a(n) = Sum_{k=0..floor(n/2)} (n+3*k)!/(k!^5 * (n-2*k)!).
a(n) = sum(k=0, n\2, (n+3*k)!/(k!^5 * (n-2*k)!));
for(n=0, M, print1(a(n),", "));

