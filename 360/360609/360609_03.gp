\\ a(n) = n! * Sum_{k=0..n} (3*k+1)^(k-1) * binomial(n+2*k,n-k)/k!.
a(n) = n! * sum(k=0, n, (3*k+1)^(k-1) * binomial(n+2*k,n-k)/k!);    
for(n=0, 18, print1(a(n), ", "))  
