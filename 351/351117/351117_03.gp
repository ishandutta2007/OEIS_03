my(N=10, x='x+O('x^N)); Vec(sum(k=0, N, k!*(k^k*x)^k/prod(j=1, k, 1-k^k*j*x)))
