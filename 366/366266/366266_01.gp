a(n) = sum(k=0, n, binomial(2*k+1,n-k) * binomial(3*k,k)/(2*k+1) );                    
for(n=0, 22, print1(a(n),", "))  

\\ A366221
b(n) = sum(k=0, n, binomial(2*k, n-k)*binomial(3*k, k)/(2*k+1));
c(n) = if(n==0, 1, b(n)+b(n-1));
for(n=0, 22, print1(a(n) - c(n),", "))
