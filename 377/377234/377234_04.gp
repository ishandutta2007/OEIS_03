a(n) = v=[1,12]; if(n<2, v[n+1], ((11*n+1)*a(n-1) - 10*(n-2)*a(n-2))/n );
for(n=0, 20, print1(a(n),", "))   