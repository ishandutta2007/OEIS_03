\\ a(n) = (1/2) * Sum_{k=0..floor(n/4)} binomial(2*n-6*k+2,2*k+1).
a(n) = sum(k=0, n\4, binomial(2*n-6*k+2, 2*k+1))/2;
for(n=0, 35, print1(a(n),", ")) 
