a(n) = sum(j=1, n, sum(k=1, n, gcd([j, k, n])^4));
for(n=1, 100, print1(a(n), ", "))

b(n) = sum(i=1, n, sum(j=1, n, sum(k=1, n, sum(l=1, n, gcd([i, j, k, l, n])^2))));
for(n=1, 30, print1(a(n)-b(n), ", "))
