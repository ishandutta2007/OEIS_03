\\ G.f.: Sum_{k>=1} x^(5*k-1) * Product_{j=1..k-1} (1-x^(4*k+j-1))/(1-x^j).
my(N=66, x='x+O('x^N)); concat([0, 0, 0], Vec(sum(k=1, N, x^(5*k-1)*prod(j=1, k-1, (1-x^(4*k+j-1))/(1-x^j)))))