M=23;

a_vector(n) = my(v=vector(n+1)); v[1]=1; for(i=1, n, v[i+1]=i!/6*sum(j=4, i, 1/(j-3)*v[i-j+1]/(i-j)!)); v;
a_vector(M)
