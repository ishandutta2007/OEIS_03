\\ a(n) = Sum_{k=0..floor(n/2)} binomial(2*k-1,n-2*k).
a(n) = sum(k=0, n\2, binomial(2*k-1, n-2*k));    
for(n=0, 41, print1(a(n),", ")) 

\\ a(n) = Sum_{k=0..floor(n/2)} binomial(2*k-2,n-2*k).
a375373(n) = sum(k=0, n\2, binomial(2*k-2, n-2*k));
b(n) = if(n==0, 1, a375373(n) + a375373(n-1));
for(n=0, 100, print1(a(n) - b(n),", "))   