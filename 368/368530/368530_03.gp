\\ a(n) = 8*a(n-1) - 22*a(n-2) + 28*a(n-3) - 17*a(n-4) + 4*a(n-5).
a(n) = my(v=[0,1,12,75,364]); if(n<5, v[n+1], 8*a(n-1) - 22*a(n-2) + 28*a(n-3) - 17*a(n-4) + 4*a(n-5) );
for(n=0, 15, print1(a(n), ", "));  
