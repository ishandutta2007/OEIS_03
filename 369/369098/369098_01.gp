my(N=30, x='x+O('x^N)); concat([0, 0, 0, 0], Vec(sum(k=1, N, k!*x^(k+3)/(1+x^4)^(k+1))))