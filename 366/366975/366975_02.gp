my(N=50, x='x+O('x^N)); concat([0, 0], Vec( sum(k=1, N, x^(3*k)/(1-x^k)^(k+1) )))   
