M=27;

a(n) = sumdiv(n, d, d*3^(d-1));
for(n=1, M, print1(a(n), ", "));