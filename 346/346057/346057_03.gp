my(N=40, x='x+O('x^N)); Vec(serlaplace(exp(-sum(k=1, N, sumdiv(k, d, 1/(d!*(k/d)^d))*x^k))))
