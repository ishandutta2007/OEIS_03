a369815(n) = sum(k=0, n\7, ((n-3*k)%4==0)*binomial((n-3*k)/4, k));
for(n=0, 50, print1(a369815(n), ", "));
\\ a(n) = A369815(7*n-4).
a(n) = if(n==0, 1, a369815(7*n-4));
for(n=0, 30, print1(a(n), ", "));