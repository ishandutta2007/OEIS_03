a(n) = my(v=[1,-2,4,-4,6,-4]); if(n<6, v[n+1], -2*a(n-1) + 4*a(n-3) + 6*a(n-4) + 4*a(n-5) + a(n-6) );                                                                     
for(n=0, 25, print1(a(n),", ")) 

   