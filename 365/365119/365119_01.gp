a(n, s=1, t=3) = sum(k=0, n, binomial(t*(n-k+1),k) * binomial(n+(s-1)*k-1,n-k)/(n-k+1) );
for(n=0, 25, print1(a(n), ", "))