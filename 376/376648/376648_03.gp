\\ a(n) = a(n-8) + a(n-9).
a(n) = my(v=[1, 0, 0, 0, 1, 0, 0, 0, 1]); if(n<9, v[n+1], a(n-8)+a(n-9));
for(n=0, 55, print1(a(n),", "))


