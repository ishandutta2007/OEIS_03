a(n) = v=[1,5]; if(n<2, v[n+1], ((2*n+3)*a(n-1) - 5*(n+3)*a(n-2))/n );
for(n=0, 28, print1(a(n),", "))   