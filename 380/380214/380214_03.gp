\\ a(n) = (1/e) * (-3)^n * n! * Sum_{k>=0} binomial(-2*k/3,n)/k!.
a(n) = round( 1/exp(1) * (-3)^n * n! * sum(k=0, 1000, binomial(-2*k/3, n)/k!) );
for(n=0, 20, print1(a(n), ", "))
