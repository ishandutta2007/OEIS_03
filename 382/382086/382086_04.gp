\\ E.g.f.: exp( Series_Reversion( x * (1-x) * exp(-x) ) ).
my(N=20, x='x+O('x^N)); Vec(serlaplace(exp( serreverse(x*(1-x)*exp(-x)) )) )

