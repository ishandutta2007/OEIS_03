M=21;

\\ a(n) = ((7*n+7)*a(n-1) - (7*n-28)*a(n-2) + (n-3)*a(n-3))/n for n > 2.
a(n) = my(v=[1,14,154]); if(n<3, v[n+1], ((7*n+7)*a(n-1) - (7*n-28)*a(n-2) + (n-3)*a(n-3))/n );
for(n=0, M, print1(a(n),", "))
