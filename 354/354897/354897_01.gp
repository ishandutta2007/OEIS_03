M=18;

a(n) = n!*sumdiv(n, d, d^n/(d!*(n/d)!));
for(n=1, M, print1(a(n), ", "));

for(n=1, 100, if(a(n)==(1+n^n), print1(n, ", ")));
