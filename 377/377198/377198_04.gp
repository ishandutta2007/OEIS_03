a(n) = v=[1,6,42]; if(n<3, v[n+1], ((7*n-1)*a(n-1) - (7*n-20)*a(n-2) + (n-3)*a(n-3))/n );
for(n=0, 20, print1(a(n),", "))   