M=15;

T(n, k) = if(k==0, (-1)^n*n!, sum(j=0, n, stirling(n, j, 2)*T(j, k-1)));
a(n) = T(n, n);
for(n=0, M, print1(a(n), ", "));
