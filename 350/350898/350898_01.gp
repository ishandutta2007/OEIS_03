my(N=99, x='x+O('x^N)); concat([0, 0, 0], Vec(sum(k=1, sqrtint(N\4), x^(4*k^2)/prod(j=1, k-1, 1-x^j))))
