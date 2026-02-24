\\ Expansion of (1/x) * Series_Reversion( x * ( Sum_{k=0..3} (-x)^k )^3 ).
my(N=30, x='x+O('x^N)); Vec(serreverse( x * sum(k=0, 3, (-x)^k)^3 ) / x)

\\ G.f.: (1/x) * Series_Reversion( x * ((1-x^4) / (1+x))^3 ).
my(N=30, x='x+O('x^N)); Vec(serreverse( x * ((1-x^4) / (1+x))^3 )/x)
