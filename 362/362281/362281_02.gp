M=19;

a(n) = n! * sum(k=0, n\2, n^k * binomial(n-k,k)/(n-k)!);
for(n=0, M, print1(a(n), ", ")); 