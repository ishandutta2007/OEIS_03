\\ a(n) = Sum_{k=0..floor(n/4)} binomial(floor(k/3),n-4*k).
a(n) = sum(k=0, n\4, binomial(k\3, n-4*k));
for(n=0, 88, print1(a(n),", "))


