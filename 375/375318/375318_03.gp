a(n) = my(v=[1,0,0,1,4,6,5]); if(n<7, v[n+1], a(n-3) + 4*a(n-4) + 6*a(n-5) + 4*a(n-6) + a(n-7) );                                                                     
for(n=0, 30, print1(a(n),", ")) 

   