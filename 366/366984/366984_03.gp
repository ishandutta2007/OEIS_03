my(N=60, x='x+O('x^N)); Vec( 1/(1-x) * sum(k=1, N, binomial(k+2,2) * x^k/(1-x^k) )) 
