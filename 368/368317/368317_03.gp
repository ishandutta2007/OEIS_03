\\ a(n) = 4^n + Sum_{k=1..n} binomial(n,k) * a(n-k).
a(n) = 4^n + sum(k=1, n, binomial(n, k) * a(n-k));
for(n=0, 20, print1(a(n), ", "));