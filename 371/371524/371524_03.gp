my(N=20, x='x+O('x^N)); Vec(2 * sum(k=0, N, (k/2+2)^(k-1) * x^k/(1 - (k/2+2)*x)^(k+1) ))