M=50;

a(n) = sumdiv(n, d, eulerphi(d)^3);
for(n=1, M, print1(a(n), ", "));
