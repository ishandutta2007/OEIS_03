\\ Expansion of e.g.f. (1/x) * Series_Reversion( x * exp(-x * (1 + x)^2) / (1 + x) ).
my(N=20, x='x+O('x^N)); Vec(serlaplace( serreverse( x * exp(-x * (1 + x)^2) / (1 + x) )/x ))


