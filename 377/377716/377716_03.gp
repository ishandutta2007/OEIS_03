\\ E.g.f.: 4/(1 + sqrt(5 - 4*exp(x)))^2.
my(N=20, x='x+O('x^N)); Vec(serlaplace( 4/(1 + sqrt(5 - 4*exp(x)))^2 ))

