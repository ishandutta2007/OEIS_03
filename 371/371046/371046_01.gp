M=20;

\\ a(n) = n! * Sum_{k=0..floor(n/3)} k^(n-3*k) * binomial(n-k+1,k)/( (n-k+1)*(n-3*k)! ).
a(n) = n! * sum(k=0, n\3, k^(n-3*k) * binomial(n-k+1,k)/( (n-k+1)*(n-3*k)! ));
for(n=0, M, print1(a(n)", "))