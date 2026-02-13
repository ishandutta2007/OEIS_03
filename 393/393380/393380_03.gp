\\ a(n) = (1/(n+1)) * Sum_{k>=0} (1/2)^(n+k+1) * binomial(n+k,k) * binomial(3*(n+k+1),n).
a(n) = sum(k=0, 2000, (1./2)^(n+k+1) * binomial(n+k,k) * binomial(3*(n+k+1), n)) / (n+1);
for(n=0, 19, print1(a(n), ", "));
for(n=0, 19, print1(round(a(n)), ", "));