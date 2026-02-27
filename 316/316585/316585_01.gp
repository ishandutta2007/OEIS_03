default(parisize, 120000000)

\\ G.f.: A(q) = ( 1 + Sum_{k>=1} (1 + q^(2*k)) * Product_{j=1..k} R_j(q) )/(1 - q*R_1(q)), where R_k(q) is defined by the continued fraction:
\\ R_k(q) = q^(2*k-1)/(1 - q^(4*k) - q^(4*k+2)/(1 - q^(4*k+4) - q^(4*k+6)/(1 - q^(4*k+8) - q^(4*k+10)/(1 - q^(4*k+12) - q^(4*k+14)/(1 ...))))).

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

M=50;
print(lista(M));
\\ v = lista(M);
\\ for(n=0, M, write("b316585_01.txt", n, " ", v[n+1]));