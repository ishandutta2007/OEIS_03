a(n) = sum(k=0, n, (-1)^(n-k) * binomial(n,k) * binomial(n+k-1,n-k) / (n-k+1) );                                
for(n=0, 31, print1(a(n),", "))   
