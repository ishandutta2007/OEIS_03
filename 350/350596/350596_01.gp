my(N=66, x='x+O('x^N)); Vec(sum(i=1, N, sum(j=1, N\i, sum(k=1, N\(i*j), x^(i*j*k)/((1-x^i)*(1-x^j)*(1-x^k))))))
