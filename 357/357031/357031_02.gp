my(A=1, n=20); for(i=1, n, A=exp((exp(x*A+x*O(x^n))-1)^2/2)); Vec(serlaplace(A))