a(n) = sum(i=1, n, sum(j=1, n, sum(k=1, n, ( gcd(i, n)/gcd([i, j, k, n]) )^2 )));   
for(n=1, 20, print1(a(n), ", "))



