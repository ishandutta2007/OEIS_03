\\ a(n) = n! * Sum_{k=1..n} 2^(n-k) * (-k)^(k-1) * binomial(n-1,k-1)/k!.
a(n) = n! * sum(k=1, n, 2^(n-k) * (-k)^(k-1) * binomial(n-1,k-1)/k!);
for(n=0, 20, print1(a(n), ", "))


