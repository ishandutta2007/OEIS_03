a(n) = sum(i=1, n, sum(j=1, n, sum(k=1, n, sum(l=1, n, ( n/gcd([i, j, k, l, n]) )^2 ))));
for(n=1, 30, print1(a(n), ", "))

\\ a(n) = Sum_{1 <= x_1, x_2, x_3, x_4 <= n} ( gcd(x_1, x_2, n)/gcd(x_1, x_2, x_3, x_4, n) )^4.
b(n) = sum(i=1, n, sum(j=1, n, sum(k=1, n, sum(l=1, n, ( gcd([i, j, n])/gcd([i, j, k, l, n]) )^4 ))));
for(n=1, 30, print1(a(n)-b(n), ", "))

