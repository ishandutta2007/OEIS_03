M=23;

a(n) = sum(k=0, n\3, k^n*binomial(n-2*k, k));
for(n=0, M, print1(a(n), ", "));

