\\ a(n) = ((n+2)/n^2) * ((2*n-1)*a(n-1) + 3*(n+1)*a(n-2)).
a(n) = my(v=[1,3]); if(n<2, v[n+1], ((n+2)/n^2) * ((2*n-1)*a(n-1) + 3*(n+1)*a(n-2)));
for(n=0, 20, print1(a(n), ", "))
