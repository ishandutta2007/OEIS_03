M=25;
N=53;

a(n) = (-1)^n + sum(k=0, (n-1)\3, a(k) * a(n-1-3*k)); 
for(n=0, M, print1(a(n),", "))

a_vector(n) = my(v=vector(n+1)); for(i=0, n, v[i+1]=(-1)^i+sum(j=0, (i-1)\3, v[j+1]*v[i-3*j])); v;
a_vector(N)