a(n, s=3, t=1, u=2) = sum(k=0, n, binomial(t*k+u*(n-k)+1, k)*binomial(s*k, n-k)/(t*k+u*(n-k)+1));
for(n=0, 23, print1(a(n), ", "))    