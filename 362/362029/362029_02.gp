my(N=100, x='x+O('x^N)); Vec(sum(k=1, N, moebius(k)^2*k*x^k)/(1+x))
