M=17;

\\ a(n) = 9^n * binomial(7*n/3 + 1/3,n)/(7*n+1).
a(n) = 9^n * binomial(7*n/3 + 1/3,n)/(7*n+1);
for(n=0, M, print1(a(n), ", ")) 