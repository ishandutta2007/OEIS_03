my(N=66, x='x+O('x^N)); Vec(sum(j=1, N, (1-x^j)*sum(k=1, N, (k*x^k)^j))/(1-x))
