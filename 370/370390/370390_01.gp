my(N=30, x='x+O('x^N)); concat([0, 0], Vec(sum(k=0, N, k!*x^k*(((1-x^2)/(1-x^3))^k-((1-x)/(1-x^2))^k))))