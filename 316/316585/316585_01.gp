default(parisize, 120000000)

lista(nn) = {
  my(q='q+O('q^(nn+1)), A=1, R=vector(nn+1), prodR=1, tmp=0);
  forstep(k=nn, 1, -1,
    tmp = q^(2*k-1)/(1 - q^(4*k) - q^(2*k+1)*tmp);
    R[k] = tmp;
  );
  for(k=1, nn,
    prodR *= R[k];
    if(valuation(prodR, 'q) > nn, break);
    A += (1 + q^(2*k))*prodR;
  );
  Vec(A/(1 - q*R[1]));
};

M=5000;
v = lista(M);
for(n=0, M, write("b316585_01.txt", n, " ", v[n+1]));