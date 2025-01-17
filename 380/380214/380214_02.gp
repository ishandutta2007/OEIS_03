M=18;

Bell(n) = n!*polcoeff(exp(exp(x+x*O(x^n))-1), n);
\\ a(n) = Sum_{k=0..n} 2^k * 3^(n-k) * |Stirling1(n,k)| * Bell(k).
a(n) = sum(k=0, n, 2^k * 3^(n-k) * abs(stirling(n, k, 1)) * Bell(k));
for(n=0, M, print1(a(n), ", "))


