\\ E.g.f.: exp( Series_Reversion( x*exp(-x)/(1+x)^2 ) ).
my(N=20, x='x+O('x^N)); Vec(serlaplace(exp( serreverse(x*exp(-x)/(1+x)^2) )) )

