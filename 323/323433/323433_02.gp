my(N=66, x='x+O('x^N)); Vec(1+sum(i=1, N, sum(j=1, N\i, x^(i*j)/prod(k=1, i*j, 1-x^k))))
