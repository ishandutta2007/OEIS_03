\\ G.f.: Sum_{k>=1} x^k * Product_{j=1..k} (1-x^(5*k+j-1))/(1-x^j).
my(N=66, x='x+O('x^N)); Vec(sum(k=1, N, x^k*prod(j=1, k, (1-x^(5*k+j-1))/(1-x^j))))
