my(N=66, x='x+O('x^N)); concat([0, 0], Vec(sum(k=1, sqrtint(N\3), x^(3*k^2)/prod(j=1, 3*k-1, 1-x^j))))
