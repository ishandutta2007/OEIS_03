my(N=40, x='x+O('x^N)); concat(0, Vec(sum(k=1, N, isprime(k)*(k*x)^k/(1-x^k)))) 
