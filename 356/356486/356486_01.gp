M=18;

a(n) = (n-1)!*sumdiv(n, d, d^n/(d-1)!);
for(n=1, M, print1(a(n), ", "));