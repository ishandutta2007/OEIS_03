\\ G.f.: (1-x-x^2) * ((1-x-x^2)^2 + 6*x^3) / ((1-x-x^2)^2 - 4*x^3)^(7/2).
my(M=30, x='x+O('x^M)); Vec( (1-x-x^2) * ((1-x-x^2)^2 + 6*x^3) / ((1-x-x^2)^2 - 4*x^3)^(7/2) )

a089627(n, k) = n!/((n-2*k)!*k!^2);
my(N=3, M=30, x='x+O('x^M), X=1-x-x^2, Y=3); Vec(sum(k=0, N\2, a089627(N, k)*X^(N-2*k)*x^(Y*k))/(X^2-4*x^Y)^(N+1/2))
