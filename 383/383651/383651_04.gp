M=24;

\\ a(n) = 3*a(n-1) + 22*a(n-2) - 24*a(n-3).
a(n) = my(v=[1, 3, 31]); if(n<3, v[n+1], 3*a(n-1) + 22*a(n-2) - 24*a(n-3));
for(n=0, M, print1(a(n), ", "));
