\\ a(n) = ((n+5)/(n*(n+4))) * ((2*n+3)*a(n-1) + 3*(n+4)*a(n-2)).
a(n) = my(v=[1,6]); if(n<2, v[n+1], ((n+5)/(n*(n+4))) * ((2*n+3)*a(n-1) + 3*(n+4)*a(n-2)));
for(n=0, 24, print1(a(n),", "))




 