M=15;

a(n) = if(n==0, 1, sum(k=1, n, k!*k^k*binomial(n-1, k-1)));
for(n=0, M, print1(a(n), ", "));