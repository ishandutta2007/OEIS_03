\\ Expansion of e.g.f. exp(x/(1 - 3*x)^3).
my(N=30, x='x+O('x^N)); Vec(serlaplace( exp(x/(1 - 3*x)^3) ))


