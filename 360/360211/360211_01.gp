M=26;

a(n) = sum(k=0, n\2, (-1)^k*binomial(2*n-3*k, n-2*k));
for(n=0, M, print1(a(n),", "))
