my(N=30, x='x+O('x^N)); Vec(serlaplace( exp(x * sum(k=0, N, binomial(4*k,k)/(3*k+1) * (4*x)^k)^4 ) ))