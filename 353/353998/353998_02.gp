M=22;

a_vector(n) = my(v=vector(n+1)); v[1]=1; for(i=1, n, v[i+1]=i!/2*sum(j=3, i, 1/(j-2)!*v[i-j+1]/(i-j)!)); v;
a_vector(M)
