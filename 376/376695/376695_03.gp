\\ a(n) = a(n-1) + a(n-7).
a(n) = my(v=[1, 1, 2, 2, 3, 3, 3]); if(n<7, v[n+1], a(n-1) + a(n-7));
for(n=0, 45, print1(a(n),", "))


