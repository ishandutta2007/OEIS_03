a(n) = if(n==1, 1, -a(n-1) + abs(moebius(n)));

for(n=1, 50, print1(a(n),", "))
