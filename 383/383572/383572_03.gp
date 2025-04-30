M=48;

\\ Diagonal of the rational function 1 / (1 - x^5 - y^5 - x^4*y^4).
a(n) = my(x='x+O('x^(n+1)), y='y+O('y^(n+1)) ); polcoef(polcoef(1 / (1 - x^5 - y^5 - x^4*y^4), n), n);
for(n=0, M, print1(a(n),", "));


