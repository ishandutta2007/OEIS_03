\\ a(n) = n! * Sum_{k=0..floor(n/2)} binomial(n+2+k,n-2*k)/k!.
a(n) = n! * sum(k=0, n\2, binomial(n+2+k, n-2*k)/k!);                                                                                 
for(n=0, 19, print1(a(n),", "))   