\\ Expansion of e.g.f. -log(1 - 2*x)/(2 * (1 + x)).
my(N=30, x='x+O('x^N)); concat(0, Vec(serlaplace( -log(1 - 2*x)/(2 * (1 + x)) )))



