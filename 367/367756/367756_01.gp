M=15;
N=21;

\\ a(0) = 1; a(n) = (n-1)! * Sum_{k=0..n-1} (k+1) * a(floor(k/4)) * a(n-1-k) / (24^floor(k/4) * floor(k/4)! * (n-1-k)!).
a(n) = if(n==0, 1, (n-1)! * sum(k=0, n-1, (k+1) * a(k\4) * a(n-1-k) / (24^(k\4) * (k\4)! * (n-1-k)!)));
for(n=0, M, print1(a(n),", "))

a_vector(n) = my(v=vector(n+1)); v[1]=1; for(i=1, n, v[i+1]=(i-1)!*sum(j=0, i-1, (j+1)*v[j\4+1]*v[i-j]/(24^(j\4)*(j\4)!*(i-1-j)!))); v;
a_vector(N)