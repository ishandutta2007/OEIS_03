a(n) = if(n==0, 1, sum(k=0, n-1, binomial(n-1,k)*d(k)));
b(n) = if(n==0, 0, sum(k=0, n-1, binomial(n-1,k)*a(k)));
c(n) = if(n==0, 0, sum(k=0, n-1, binomial(n-1,k)*b(k)));
d(n) = if(n==0, 0, sum(k=0, n-1, binomial(n-1,k)*c(k)));

for(n=0, 20, print1(b(n),", "))      


