my(N=20, x='x+O('x^N)); concat(0, Vec(serlaplace(sum(k=1, N, (1 - exp(-k*x))^k / k^3 ))))