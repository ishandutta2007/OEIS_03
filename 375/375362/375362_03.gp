a(n) = my(v=[1,0,4,6,20]); if(n<5, v[n+1], 4*a(n-2) + 6*a(n-3) + 4*a(n-4) + a(n-5) );                                                                     
for(n=0, 30, print1(a(n),", ")) 

   