M=18;

a(n)=local(A=1); for(i=1, n, A=1+2*(1-exp(-x+x*O(x^n)))*subst( A, x, 1-exp(-x+x*O(x^n)) )); n!*polcoeff(A, n);
for(n=0, M, print1(a(n), ", "));