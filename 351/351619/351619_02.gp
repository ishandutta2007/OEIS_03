my(N=99, x='x+O('x^N)); concat(0, Vec(sum(k=1, N, isprime(k)*(-x)^k/(1-x^k))))
