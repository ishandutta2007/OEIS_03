\\ a(n) = n! * Sum_{k=0..floor(n/2)} (-1)^k * binomial(n/2-k,k)/(n-2*k)!.
a(n) = n! * sum(k=0, n\2, (-1)^k * binomial(n/2-k,k)/(n-2*k)!);
for(n=0, 23, print1(a(n),", "))
