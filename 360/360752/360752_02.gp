M=41;

a(n) = sumdiv(n, d, 2^(n-d)*binomial(d, n/d-1));
for(n=1, M, print1(a(n), ", "));
for(n=1, M, print1(a(n)-1, ", "));