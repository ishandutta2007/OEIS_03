a(n) = n! * sum(k=0, n, (n+1)^(k-1) * binomial(n+k-1,n-k)/k! );
for(n=0, 17, print1(a(n),", "))      
 
