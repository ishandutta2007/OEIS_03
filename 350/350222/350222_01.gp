M=43;

a(n) = sum(k=1, n, (-1)^(k+1)*(n^3\k^3));
for(n=1, M, print1(a(n), ", "));
