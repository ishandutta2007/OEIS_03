\\ Sum_{k in Z} x^k / (1 - x^(7*k+5)).
my(N=50, x='x+O('x^N)); Vec(sum(k=-N, N, x^(k) / (1 - x^(7*k+5))))

\\ G.f.: Sum_{k in Z} x^(5*k) / (1 - x^(7*k+1)).
my(N=50, x='x+O('x^N)); Vec(sum(k=-N, N, x^(5*k) / (1 - x^(7*k+1))))
