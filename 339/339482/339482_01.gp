M=27;

a(n) = sumdiv(n, d, d^(n-d+1)*binomial(d+n/d-2, d-1));
for(n=1, M, print1(a(n), ", "));
