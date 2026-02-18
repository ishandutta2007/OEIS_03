\\ G.f.: Sum_{k>=0} x^(7*k) / ( (Product_{j=1..k} (1-x^(2*j))) * (Product_{j=1..5*k} (1-x^(2*j))) ).
my(N=70, x='x+O('x^N)); Vec(sum(k=0, N, x^(7*k)/(prod(j=1, k, (1-x^(2*j)))*prod(j=1, 5*k, (1-x^(2*j))))))
my(N=70, x='x+O('x^N)); sum(k=0, N, x^(7*k)/(prod(j=1, k, (1-x^(2*j)))*prod(j=1, 5*k, (1-x^(2*j)))))