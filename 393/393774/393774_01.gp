lista(nn) = {
  my(q='q+O('q^(nn+1)), R=1);
  forstep(k=nn, 0, -1, R = 1 - q^(2*k+2) - q^(2*k+3)/R);
  Vec(1/R, -nn-1);
};
print(lista(46));