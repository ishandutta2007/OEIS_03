M=65;

a(n) = sumdiv(n, d, (d%2)*binomial(d+n/d-1, d));
for(n=1, M, print1(a(n),", "));
