M=10;

T(n, k) = round(((k+sqrt(k))^n+(k-sqrt(k))^n))/2;
for(n=0, M, for(k=0, n, print1(T(k, n-k), ", ")));


