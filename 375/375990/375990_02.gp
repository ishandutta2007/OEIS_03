\\ Expansion of e.g.f. (1 + 2 * log(1 - x))^2.
my(N=25, x='x+O('x^N)); Vec(serlaplace((1 + 2 * log(1 - x))^2))


