a(n, r=3, s=1, t=1, u=0) = r*n!*sum(k=0, n, (t*k+u*(n-k)+r)^(k-1)*binomial(s*k, n-k)/k!);
for(n=0, 17, print1(a(n), ", "))

