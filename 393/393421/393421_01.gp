\\ Expansion of (1/x) * Series_Reversion( x * (5 * (1 - x)^2 - 4) ).
my(N=30, x='x+O('x^N)); Vec(serreverse(x*(5*(1-x)^2-4))/x)