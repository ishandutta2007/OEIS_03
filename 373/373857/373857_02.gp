my(N=20, x='x+O('x^N)); concat(0, Vec(serlaplace(sum(k=1, N, log(1 + k*x)^k / k ))))