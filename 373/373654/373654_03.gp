a(n) = my(v=[1, 1, 1, 1, 3, 5, 7]); if(n<7, v[n+1], a(n-1) + 2*a(n-3) - a(n-6) );
for(n=0, 25, print1(a(n),", "))