my(N=30, x='x+O('x^N)); Vec(serlaplace(sum(k=1, N, (1-x^k)*(exp(x^k)-1)/k!)/(1-x)))