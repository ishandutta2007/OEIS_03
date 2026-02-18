\\ G.f.: Sum_{k>=0} x^(3*k) / ( (Product_{j=1..k} (1-x^j)) * (Product_{j=1..4*k} (1-x^j)) ).
my(N=70, x='x+O('x^N)); Vec(sum(k=0, N, x^(3*k)/(prod(j=1, k, (1-x^j))*prod(j=1, 4*k, (1-x^j)))))
my(N=70, x='x+O('x^N)); sum(k=0, N, x^(3*k)/(prod(j=1, k, (1-x^j))*prod(j=1, 4*k, (1-x^j))))