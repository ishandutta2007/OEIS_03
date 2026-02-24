\\ Expansion of (1/x) * Series_Reversion( x * ( Sum_{k=0..5} (-x)^k )^2 ).
my(N=30, x='x+O('x^N)); Vec(serreverse( x * sum(k=0, 5, (-x)^k)^2 ) / x)

\\ G.f.: (1/x) * Series_Reversion( x * ((1-x^6) / (1+x))^2 ).
my(N=30, x='x+O('x^N)); Vec(serreverse( x * ((1-x^6) / (1+x))^2 )/x)
