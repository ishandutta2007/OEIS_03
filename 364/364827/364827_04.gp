a(n) = if(n==0, 1, sum(k=1, n, 2^k * 3^(n-k) * binomial(n,k) * binomial(5*n,k-1) )/n); 
for(n=0, 17, print1(a(n),", "))  

