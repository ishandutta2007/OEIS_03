M=21;

\\ a(n) = ((7*n+3)*a(n-1) - (7*n-24)*a(n-2) + (n-3)*a(n-3))/n for n > 2.
a(n) = my(v=[1,10,90]); if(n<3, v[n+1], ((7*n+3)*a(n-1) - (7*n-24)*a(n-2) + (n-3)*a(n-3))/n );
for(n=0, M, print1(a(n),", "))
