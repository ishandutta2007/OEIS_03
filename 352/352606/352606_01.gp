my(N=40, x='x+O('x^N)); v=Vec(serlaplace(tan(cos(x)-1))); concat(0, vector(#v\2, k, v[2*k-1]))

