a(n) = v=[1,15]; if(n<2, v[n+1], ((11*n+4)*a(n-1) - 10*(n-2)*a(n-2))/n );
for(n=0, 20, print1(a(n),", "))   