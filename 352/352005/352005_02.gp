my(N=40, x='x+O('x^N)); Vec(serlaplace(exp(sum(k=1, N, sumdiv(k, d, isprime(d)*(-1)^(k/d+1)*(k-1)!/(d-1)!)*x^k/k!))))
