a(n) =my(v=[1, 0, 1]); if(n<3, v[n+1], (n-1)/2 * (2*a(n-2) + 3*(n-2)*a(n-3)) );  
for(n=0, 24, print1(a(n),", "))  
