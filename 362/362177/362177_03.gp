a(n) = n! * sum(k=0, n\2, (-3)^k * binomial(n-k,k)/(n-k)!);

for(n=0, 22, print1(a(n),", "))          