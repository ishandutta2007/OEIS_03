a(n) = v=[1,2,-2,-8,6]; if(n<5,v[n+1], -( (n-3)*a(n-1) + (6*n-6)*a(n-2) + 10*(n-3)*a(n-3) + 5*(n-4)*a(n-4) + (n-5)*a(n-5) )/n);
for(n=0, 20, print1(a(n),", "))  
