a(n) = if(n==0, 1, sum(k=0, (n-4)\5, binomial(n,5*k+4) * a(n-5*k-4) ));                   
for(n=0, 29, print1(a(n),", "))  
