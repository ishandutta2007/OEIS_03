\\ G.f.: Sum_{k>=0} x^(4*k) / ( (Product_{j=1..2*k} (1-x^j)) * (Product_{j=1..3*k} (1-x^j)) ).
my(N=70, x='x+O('x^N)); Vec(sum(k=0, N, x^(4*k)/(prod(j=1, 2*k, (1-x^j))*prod(j=1, 3*k, (1-x^j)))))
my(N=70, x='x+O('x^N)); sum(k=0, N, x^(4*k)/(prod(j=1, 2*k, (1-x^j))*prod(j=1, 3*k, (1-x^j))))