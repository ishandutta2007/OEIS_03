\\ Expansion of (1/x) * Series_Reversion( x * (4 * (1 - x)^2 - 3) ).
my(N=30, x='x+O('x^N)); Vec(serreverse(x*(4*(1-x)^2-3))/x)