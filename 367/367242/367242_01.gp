a(n, s=2, t=0, u=3) = sum(k=0, n, binomial(t*k+u*(n-k)+1, k)*binomial(n+(s-1)*k-1, n-k)/(t*k+u*(n-k)+1));
for(n=0, 24, print1(a(n), ", "))    