\\ E.g.f.: (1/x) * Series_Reversion( x/(1 - x*log(1 - x)) ). 
my(N=25, x='x+O('x^N)); Vec(serlaplace( serreverse( x/(1 - x*log(1 - x)) )/x ))

