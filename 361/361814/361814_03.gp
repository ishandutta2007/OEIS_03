a(n) = v=[1,2,16,100,660,4482]; if(n<6,v[n+1], 2 * ( (2*n-1)*a(n-1) + 5*(2*n-2)*a(n-2) + 10*(2*n-3)*a(n-3) + 10*(2*n-4)*a(n-4) + 5*(2*n-5)*a(n-5) + (2*n-6)*a(n-6) )/n);           
for(n=0, 20, print1(a(n),", ")) 