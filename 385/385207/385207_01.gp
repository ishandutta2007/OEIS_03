M=18;

\\ a(n) = 49^n * binomial(6*n/7+1/7,n)/(6*n+1).
a(n) = 49^n * binomial(6*n/7+1/7,n)/(6*n+1);
for(n=0, M, print1(a(n), ", ")) 

\\ a(7*n+1) = 0 for n > 0.
for(n=0, 100, if(a(7*n+1) != 0, print1(n, ", ")));