a(n) = my(v=[1,0,1,0,1,0,1,1]); if(n<7, v[n+1], a(n-2) + a(n-7) );
for(n=0, 35, print1(a(n), ", "));