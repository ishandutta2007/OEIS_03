a(n) = polcoef(1/sqrt(1 - 4*x*(1+x)^n + x*O(x^n)), n);
for(n=0, 19, print1(a(n),", ")) 

