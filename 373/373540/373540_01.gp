\\ a(n) = n! * Sum_{k=0..floor(n/4)} (-1)^k * binomial(n/4-1,k)/(n-4*k)!.
a(n) = n! * sum(k=0, n\4, (-1)^k * binomial(n/4-1,k)/(n-4*k)!);
for(n=0, 22, print1(a(n),", "))
