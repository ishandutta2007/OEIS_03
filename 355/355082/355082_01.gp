M=12;

a_vector(n) = my(v=vector(n+1)); v[1]=1; for(i=1, n, v[i+1]=sum(j=0, i-1, 4^j*binomial(i-1, j)*v[j+1])); v;
a_vector(M)