\\ a(n) = 5*a(n-1) - 4*a(n-2) + a(n-3) - 5*a(n-4) + 4*a(n-5).
a(n) = my(v=[0,0,0,1,5]); if(n<5, v[n+1], 5*a(n-1) - 4*a(n-2) + a(n-3) - 5*a(n-4) + 4*a(n-5));
for(n=0, 20, print1(a(n), ", "));


 