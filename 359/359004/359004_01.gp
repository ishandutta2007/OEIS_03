M=43;

a(n) = sumdiv(n, d, d^(n/d-1)*(n/d)^(d-1));
for(n=1, M, print1(a(n), ", "));