M=18;

a_vector(n) = my(v=vector(n+1)); v[1]=0; v[2]=1; for(i=2, n, v[i+1]=(3*i-2)*v[i]-3*(i-1)*v[i-1]+1); v;
a_vector(M)
