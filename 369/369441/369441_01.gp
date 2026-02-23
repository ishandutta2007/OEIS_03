my(N=30, x='x+O('x^N)); Vec(serreverse( x / ((1+x)^2 * (1+x^2)^2) )/x)

\\ G.f.: B(x)^2, where B(x) is the g.f. of A186241.
my(N=30, x='x+O('x^N)); Vec((serreverse( x / ((1+x)^2 * (1+x^2)^2) )/x)^(1/2))
