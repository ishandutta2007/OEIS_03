M=20;

a(n) = sum(k=0, n\2, k^n*binomial(n-k-1, n-2*k));
for(n=0, M, print1(a(n), ", "));

