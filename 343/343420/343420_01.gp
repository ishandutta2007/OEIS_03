M=26;

a(n) = my(A=1+O(x)); for(i=1, n, A=1-((n-i+1)*x)^(n-i+1)/A); polcoef(1/A, n);
for(n=0, M, print1(a(n), ", "));
