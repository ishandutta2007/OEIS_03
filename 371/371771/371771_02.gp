\\ a(n) = [x^n] 1/((1-x^3) * (1-x)^(3*n)).
a(n) = polcoef(1/((1-x^3) * (1-x)^(3*n) + x*O(x^n)), n);
for(n=0, 19, print1(a(n), ", "))