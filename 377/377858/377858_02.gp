\\ G.f.: x * (9 + 45*x + 11*x^2 - x^3)/(1 - x)^5.
my(N=40, x='x+O('x^N)); concat(0, Vec( x * (9 + 45*x + 11*x^2 - x^3)/(1 - x)^5 ))


